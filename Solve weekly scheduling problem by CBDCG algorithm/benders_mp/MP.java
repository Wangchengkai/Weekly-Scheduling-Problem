package com.ben.benders_mp;

import com.ben.enums.LocalSearchSolStatus;
import com.ben.benders_mp.util.Coder;
import com.ben.params;
import com.ben.enums.SolStatus;
import gurobi.*;
import lombok.Data;
import lombok.extern.slf4j.Slf4j;

import java.io.*;
import java.util.*;

@Slf4j
@Data
public class MP {
    private GRBModel model;
    private GRBEnv envForLocalSearch;               // the GRB environment for local search
    private GRBVar[] twd;                           // GRB variable: the total working duration of each day
    private GRBVar[] costOfAddedPhysicianPerDay;    // the cost of adding seconded physicians per day
    private double[] costOfAddedPhysicianValPerDay; // solution result: the cost of adding seconded physicians per day
    private double[] twdVal;                        // solution result: the total working duration of each day
    private GRBVar[][] p_dt;                        // GRB variable: the number of physicians
    private GRBVar[][] r;                           // GRB variable: the selection of blocks
    private GRBVar[][][] p_dtk;                     // GRB variable: the indicator variable of p_dt
    private GRBVar[][] L_dt;                        // GRB variable: the queue length of each day and each period
    private double[][] range;                       // the upper bound of service stress rho
    private double[][] slope;                       // the slope of piecewise linear function
    private double[][] inter1;                      // the upper bound of service stress rho
    private double[][] slope1;                      // the slope of piecewise linear function
    private GRBVar[][][] z;                         // the indicator variable of added cut
    private GRBVar[][] zz;

    private double objVal = 0;                      // the objective value
    private double[][] p;                           // the number of physicians in each day and each period
    private double[][] queueLengh;                  // queue length in each day and each period
    private int genCutCount;                        // the number of generated cuts
    private GRBGenConstr[] latestGenConstrs;        // the latest generated constraints
    private List<int[]> shifts;                     // all possible shifts
    private Set<String> addedCuts;                  // the set of added cuts, to avoid adding duplicate cuts
    private Map<String, LocalSearchSolStatus> searchedCuts;// the set of searched cuts, to avoid evaluating duplicate cuts
    private Set<String> addedCutsWeekly;            // the set of added cuts(weekly), to avoid adding duplicate cuts
    private Map<String, LocalSearchSolStatus> searchedCutsWeekly;// the set of searched cuts(weekly), to avoid evaluating duplicate cuts
    private Set<String> addedCutsForPhy;            // the set of added cuts, to avoid adding duplicate cuts
    private Map<String, LocalSearchSolStatus> searchedCutsWeeklyForPhy;// the set of searched cuts, to avoid evaluating duplicate cuts

    private GRBVar[][] cost_dt;
    private GRBVar[] cost_d;
    private double[] cost_d_Val;

    public MP(GRBEnv env) {
        try {
            env.set(GRB.DoubleParam.MIPGap, 0.01);
            envForLocalSearch = new GRBEnv();
            model = new GRBModel(env);
            model.set(GRB.IntParam.OutputFlag, 0);// close the log
            addedCuts = new HashSet<>();
            addedCutsWeekly = new HashSet<>();
            addedCutsForPhy =new HashSet<>();

            searchedCuts = new HashMap<>();
            p = new double[params.D][params.T];
            genCutCount = 0;
            readParams(params.PWL_PARAMS_PATH);// read the parameters of piecewise linear function
            readKandB(params.PWL_PARAMS_PATH_wck);
            initModel();
        } catch (GRBException | IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * read the parameters of piecewise linear function
     * @param filename  the file of parameters
     */
    public void readParams(String filename) throws IOException {
        File file = new File(filename);
        range = new double[params.N][params.PARTITIONS];
        slope = new double[params.N][params.PARTITIONS];
        BufferedReader br = new BufferedReader(new FileReader(file));
        for (int k = 0; k < params.N; k++) {
            for (int v = 0; v < params.PARTITIONS; v++) {
                String[] tmp = br.readLine().split(" ");
                range[k][v] = Double.parseDouble(tmp[0]);
                slope[k][v] = Double.parseDouble(tmp[1]) * 10000;
            }
            br.readLine();
        }
    }

    public void readKandB(String filename) throws IOException{
        File file = new File(filename);
        slope1 = new double[params.N][params.PARTITIONS];
        inter1 = new double[params.N][params.PARTITIONS];
        BufferedReader br = new BufferedReader(new FileReader(file));
        for (int k = 0; k < params.N; k++) {
            for (int v = 0; v < params.PARTITIONS; v++) {
                String[] tmp = br.readLine().split(" ");
                slope1[k][v] = Double.parseDouble(tmp[0]);
                inter1[k][v] = Double.parseDouble(tmp[1]);
            }
            br.readLine();
        }
    }

    // the call function of piecewise linear function
    public double calQueueLengthPWL(double rho, int s) {
        double ql = 0;
        for (int j = 0; j < params.PARTITIONS; j++) {
            if (range[s - 1][j] <= rho) {
                ql += range[s - 1][j] * slope[s - 1][j];
                rho -= range[s - 1][j];
            } else {
                ql += rho * slope[s - 1][j];
                break;
            }
        }
        return ql;
    }

    /**
     * init the model of the master problem: variables + constraints
     * name format:
     * var: x[i]
     * constr: c(i)
     */
    private void initModel() {
        genAllShifts();
        try {
            // add vars
            costOfAddedPhysicianPerDay=new GRBVar[params.D];
            for (int d = 0; d < params.D; d++) {
                costOfAddedPhysicianPerDay[d] = model.addVar(0,  params.N*params.PuhishCostForAddedPhysician, 1, GRB.CONTINUOUS, String.format("costForPhyDaily[%d]", d + 1));
            }

            twd = new GRBVar[params.D];//作为一个一维变量的定义和赋值。
            for (int d = 0; d < params.D; d++) {
                twd[d] = model.addVar(params.T, 5 * params.T, 1, GRB.INTEGER, String.format("twd[%d]", d + 1));
            }

            cost_dt=new GRBVar[params.D][params.T];
            for (int d = 0; d < params.D; d++)
                for (int t = 0; t < params.T; t++)
                    cost_dt[d][t]= model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, String.format("cost_dt[%d,%d]", d + 1, t + 1));

            cost_d=new GRBVar[params.D];
            for (int d = 0; d < params.D; d++) {
                cost_d[d]= model.addVar( 0, GRB.INFINITY, 1, GRB.CONTINUOUS, String.format("cost_d[%d]", d + 1));
            }

            p_dt = new GRBVar[params.D][params.T];   // the number of physicians
            for (int d = 0; d < params.D; d++)
                for (int t = 0; t < params.T; t++)
                    p_dt[d][t] = model.addVar(1, params.N_2 - 1, 0, GRB.INTEGER, String.format("p[%d,%d]", d + 1, t + 1));


            p_dtk = new GRBVar[params.D][params.T][params.N_2];
            for (int d = 0; d < params.D; d++)
                for (int t = 0; t < params.T; t++)
                    for (int k = 0; k < params.N_2; k++)
                        p_dtk[d][t][k] = model.addVar(0, 1, 0, GRB.INTEGER, String.format("p[%d,%d,%d]", d + 1, t + 1, k + 1));

            L_dt = new GRBVar[params.D][params.T]; // queue length
            for (int d = 0; d < params.D; d++)
                for (int t = 0; t < params.T; t++)
                    L_dt[d][t] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, String.format("L[%d,%d]", d + 1, t + 1));


            // change the definition of the variable
            GRBVar[][] rho1 = new GRBVar[params.D][params.T];  //  the self variable of piecewise linear function
            for (int d = 0; d < params.D; d++)
                for (int t = 0; t < params.T; t++){
                    rho1[d][t] = model.addVar(0, 1, 0, GRB.CONTINUOUS, String.format("rho[%d,%d]", d + 1, t + 1));
                }

            // about the schedule constraints
            r = new GRBVar[params.D][shifts.size()];
            for (int d = 0; d < params.D; d++)
                for (int j = 0; j < shifts.size(); j++)
                    r[d][j] = model.addVar(0, GRB.INFINITY, 0, GRB.INTEGER, String.format("r[%d,%d]", d + 1, j + 1));

            //about the new added cut
            z = new GRBVar[params.D][params.T][params.MaxIteration];
            for (int d = 0; d < params.D; d++)
                for (int t = 0; t < params.T; t++)
                    for(int iter=0;iter<params.MaxIteration;iter++)
                        z[d][t][iter] = model.addVar(0, 1, 0, GRB.INTEGER, String.format("z[%d,%d]", d + 1, t + 1));

            zz = new GRBVar[params.T][params.D * params.MaxIteration];
            for (int t = 0; t < params.T; t++)
                for(int iter=0;iter<params.D * params.MaxIteration;iter++)
                    zz[t][iter] = model.addVar(0, 1, 0, GRB.INTEGER, String.format("zz[%d]", t + 1));

            // add constraints
            model.set(GRB.IntAttr.ModelSense, GRB.MINIMIZE);  //目标函数在p_t变量的定义中已经设置好了

            // total working time per day
            for (int d = 0; d < params.D; d++) {
                GRBLinExpr expr = new GRBLinExpr();
                expr.addTerm(1, twd[d]);
                for (int t = 0; t < params.T; t++) {
                    expr.addTerm(-1, p_dt[d][t]);
                }
                model.addConstr(expr, GRB.GREATER_EQUAL, 0, String.format("twd(%d)", d + 1));
            }

            // bind p and schedule
            for (int d = 0; d < params.D; d++) {
                for (int t = 0; t < params.T; t++) {
                    GRBLinExpr bindShift = new GRBLinExpr();
                    for (int j = 0; j < shifts.size(); j++) {
                        double n = shifts.get(j)[0];
                        double endT = shifts.get(j)[1];
                        if (t >= endT - n + 1 && t <= endT)
                            bindShift.addTerm(1.0, r[d][j]);
                    }
                    model.addConstr(bindShift, GRB.EQUAL, p_dt[d][t], "bindShift(" + (d + 1) + "," + (t + 1) + ")");
                }
            }

            // one valid inequality
            for (int d = 0; d < params.D; d++){
                GRBLinExpr AboutMinPhyUsed=new GRBLinExpr();
                for (int j = 0; j < shifts.size(); j++){
                    AboutMinPhyUsed.addTerm(0.5, r[d][j]);
                }
                for(int m=params.IndexOfFinicPhysician;m<params.N;m++){
                    AboutMinPhyUsed.addTerm(-1/params.PuhishCostForAddedPhysician, costOfAddedPhysicianPerDay[d]);
                }
                model.addConstr(AboutMinPhyUsed, GRB.LESS_EQUAL, params.IndexOfFinicPhysician, "VI"+(d+1));

            }

            for (int d = 0; d < params.D; d++) {
                for (int t = 0; t < params.T; t++) {     // bind queue length and the indicator variable
                    GRBLinExpr pIndBind = new GRBLinExpr();
                    for (int k = 0; k < params.N_2; k++) {
                        pIndBind.addTerm(k + 1, p_dtk[d][t][k]);
                    }
                    model.addConstr(pIndBind, GRB.EQUAL, p_dt[d][t], "pIndBind(" + (d + 1) + "," + (t + 1) + ")");
                }
            }

            for (int d = 0; d < params.D; d++) {
                for (int t = 0; t < params.T; t++) {     // the constraint of the indicator variable of the queue length
                    GRBLinExpr pInd = new GRBLinExpr();
                    for (int k = 0; k < params.N_2; k++) {
                        pInd.addTerm(1, p_dtk[d][t][k]);
                    }
                    model.addConstr(pInd, GRB.EQUAL, 1, "pInd(" + (d + 1) + "," + (t + 1) + ")");
                }
            }

            // the constraint of the maximum queue length
            if(params.isHardContrOfQueueLength){
                for (int d = 0; d < params.D; d++) {
                    for (int t = 0; t < params.T; t++) {
                        GRBLinExpr maxQL = new GRBLinExpr();
                        maxQL.addTerm(1.0, L_dt[d][t]);
                        model.addConstr(maxQL, GRB.LESS_EQUAL, params.Lr, "maxQL(" + (d + 1) + "," + (t + 1) + ")");
                    }
                }
            }

            // the constraint of staffing
            for (int d = 0; d < params.D; d++) {
                for (int t = 0; t < params.T-params.LEN_NIGHT_SCH; t++) {
                    model.addConstr(p_dt[d][t], GRB.GREATER_EQUAL, params.minOfStaffig, "");
                }
            }

            // the punish cost
            for (int d = 0; d < params.D; d++) {
                for (int t = 0; t < params.T; t++) {
                    GRBLinExpr costInOneHour = new GRBLinExpr();
                    costInOneHour.addTerm(1.0, cost_dt[d][t]);
                    costInOneHour.addTerm(- params.PunishCost, L_dt[d][t]);
                    model.addConstr(costInOneHour, GRB.GREATER_EQUAL, - params.Lr * params.PunishCost, "pubish(" + (d + 1)+","+(t+1) + ")");
                }
            }
            for (int d = 0; d < params.D; d++){
                GRBLinExpr costInOneDay = new GRBLinExpr();
                costInOneDay.addTerm(1, cost_d[d]);
                for (int t = 0; t < params.T; t++) {
                    costInOneDay.addTerm(-1, cost_dt[d][t]);
                }
                model.addConstr(costInOneDay, GRB.EQUAL, 0, "pubish(" + (d + 1) + ")");
            }

            for (int d = 0; d < params.D; d++) {
                for (int t = 0; t < params.T; t++) {     // compute the queue length
                    for (int k = 0; k < params.N_2; k++) {
                        for (int v = 0; v < params.PARTITIONS; v++) {
                            GRBLinExpr QL_left = new GRBLinExpr();
                            QL_left.addTerm(slope1[k][v], rho1[d][t]);
                            QL_left.addTerm(slope1[k][v], p_dtk[d][t][k]);
                            QL_left.addTerm(inter1[k][v], p_dtk[d][t][k]);
                            QL_left.addTerm(-1.0, L_dt[d][t]);
                            model.addConstr(QL_left, GRB.LESS_EQUAL, slope1[k][v],"QL_left("+(d + 1) + "," + (t + 1) + "," + (k + 1) + "," + (v + 1) + ")");
                        }
                    }
                }
            }

            // Conservation Of Flow
            for (int d = 0; d < params.D; d++) {

                // t=1
                GRBLinExpr COF = new GRBLinExpr();
                for (int k = 0; k < params.N_2; k++){
                    COF = new GRBLinExpr();
                    COF.addTerm(1.0, L_dt[d][0]);
                    COF.addTerm((k+1)*params.MU * params.DELTA, rho1[d][0]);
                    COF.addTerm(-params.M, p_dtk[d][0][k]);
                    model.addConstr(COF, GRB.GREATER_EQUAL, params.LAM[d][0] * params.DELTA-params.M, "COF(" + (d + 1) + "," + 1 +","+ (k+1)+ ")");
                }

                COF = new GRBLinExpr();
                for (int k = 0; k < params.N_2; k++){
                    COF = new GRBLinExpr();
                    COF.addTerm(1.0, L_dt[d][0]);
                    COF.addTerm((k+1)*params.MU * params.DELTA, rho1[d][0]);
                    COF.addTerm(params.M, p_dtk[d][0][k]);
                    model.addConstr(COF, GRB.LESS_EQUAL, params.LAM[d][0] * params.DELTA+params.M, "COF(" + (d + 1) + "," + 1 +","+ (k+1)+ ")");
                }

                // t>1
                for (int t = 1; t < params.T; t++) {
                    for (int k = 0; k < params.N_2; k++) {
                        COF = new GRBLinExpr();
                        COF.addTerm(1.0, L_dt[d][t]);
                        COF.addTerm(-1.0, L_dt[d][t - 1]);
                        COF.addTerm((k+1)*params.MU * params.DELTA, rho1[d][t]);
                        COF.addTerm(-params.M, p_dtk[d][t][k]);
                        model.addConstr(COF, GRB.GREATER_EQUAL, params.LAM[d][t] * params.DELTA-params.M, "COF(" + (d + 1) + "," + (t+1) +","+ (k+1)+ ")");
                    }
                    for (int k = 0; k < params.N_2; k++) {
                        COF = new GRBLinExpr();
                        COF.addTerm(1.0, L_dt[d][t]);
                        COF.addTerm(-1.0, L_dt[d][t - 1]);
                        COF.addTerm((k+1)*params.MU * params.DELTA, rho1[d][t]);
                        COF.addTerm(params.M, p_dtk[d][t][k]);
                        model.addConstr(COF, GRB.LESS_EQUAL, params.LAM[d][t] * params.DELTA+params.M, "COF(" + (d + 1) + "," + (t+1) +","+ (k+1)+ ")");
                    }
                }
            }

            //  the number of doctors should be the same during the night shift
            for (int d = 0; d < params.D; d++) {
                for (int t = 16; t < params.T; t++) {
                    model.addConstr(p_dt[d][t], GRB.EQUAL, p_dt[d][params.T - 1], String.format("nightEq(%d,%d)", d + 1, t + 1));
                }
            }

            //  remove the repeated solutions (a valid inequalities)
            for (int i = 0; i < shifts.size(); i++) {
                for (int j = i + 1; j < shifts.size(); j++) {
                    if (eitherOr(i, j)) {
                        for (int d = 0; d < params.D; d++) {
                            double BigM = 40;
                            GRBVar w = model.addVar(0, GRB.INFINITY, 0, GRB.INTEGER, String.format("w[%d,%d,%d]", d + 1, i, j));
                            GRBLinExpr expr = new GRBLinExpr();
                            expr.addTerm(1, r[d][i]);
                            expr.addTerm(-BigM, w);
                            model.addConstr(expr, GRB.LESS_EQUAL, 0, String.format("con1(%d,%d,%d)", d + 1, i, j));
                            GRBLinExpr expr2 = new GRBLinExpr();
                            expr2.addTerm(1, r[d][j]);
                            expr2.addTerm(BigM, w);
                            model.addConstr(expr2, GRB.LESS_EQUAL, BigM, String.format("con2(%d,%d,%d)", d + 1, i, j));
                        }
                    }
                }
            }
        } catch (GRBException e) {
            e.printStackTrace();
        }
    }

    private void genAllShifts() {
        // {n,t}:{the length of shift, the end time of shift}
        shifts = new ArrayList<>();
        // 空班次：即不用该医生
        shifts.add(new int[]{0, 0});
        // 夜班
        shifts.add(new int[]{8, 23});
        for (int n = params.S_LB; n <= params.S_UB; n++)
            for (int t = n - 1; t < params.T - 8; t++)
                shifts.add(new int[]{n, t});
    }

    // two shifts cannot exist at the same time
    private boolean eitherOr(int i, int j) {
        int[] shiftOne = shifts.get(i);
        int[] shiftTwo = shifts.get(j);
        return (shiftOne[1] - shiftOne[0] - (shiftTwo[1] - shiftTwo[0])) * (shiftOne[1] - shiftTwo[1]) < 0;
    }


    public void saveResult() throws IOException {
        BufferedWriter br = new BufferedWriter(new FileWriter("data/result.txt", true));
        br.write("*****************  MP RESULT  *****************\n");
        br.write("the objective value: " + objVal + "\n");
        br.write("the number of physician in each period: ");
        for (int d = 0; d < params.D; d++) {
            for (int t = 0; t < params.T; t++) {
                br.write((int) Math.round(p[d][t]) + ", ");
            }
            br.write("\n");
        }

        br.write("the queue length in each period \n");
        for (int d = 0; d < params.D; d++) {
            for (int t = 0; t < params.T; t++) {
                br.write(String.format("%.2f, ", queueLengh[d][t]));
            }
            br.write('\n');
        }
        br.write("\n");
        br.close();
    }

    public void addCut(double[] twd, double[] costOfAddedPhysicianPerDay) {
        for (int d = 0; d < params.D; d++) {
            addCutToModel_wck(p, d, twd[d],costOfAddedPhysicianPerDay[d]);
        }
    }

    public void addCutToModel_wck(double[][] p_dt_double, int d, double lb,double lbForPhy) {
        String cutEncoded = Coder.encodeCut(p_dt_double[d], d, lb, lbForPhy);
        if (addedCuts.contains(cutEncoded))//ifdef-notdef
            return;
        addedCuts.add(cutEncoded);

        try {
            if(addedCuts.size()%100==0)
                System.out.println("cut size is N * 100: " + addedCuts.size());
            if(addedCuts.size()>params.MaxIteration*params.D){
                System.out.println("cut size is too large: " + addedCuts.size());
            }

            GRBLinExpr cut = new GRBLinExpr();
            cut.addTerm(1, twd[d]);
            for (int t = 0; t < params.T; t++) {
                cut.addTerm(-lb, zz[t][addedCuts.size()-1]);
            }
            model.addConstr(cut, GRB.GREATER_EQUAL, (1 - params.T) * lb, "Gen(" + (genCutCount++) + ")");

            //PhyCost
            GRBLinExpr cut2 = new GRBLinExpr();
            cut2.addTerm(1, costOfAddedPhysicianPerDay[d]);
            for (int t = 0; t < params.T; t++) {
                cut2.addTerm(-lbForPhy, zz[t][addedCuts.size() - 1]);
            }
            model.addConstr(cut2, GRB.GREATER_EQUAL, (1 - params.T) * lbForPhy, "GenForPhy(" + (genCutCount++) + ")");

            for (int t = 0; t < params.T; t++) {
                GRBLinExpr cut1 = new GRBLinExpr();
                cut1.addTerm(params.M, zz[t][addedCuts.size()-1]);
                cut1.addTerm(-1, p_dt[d][t]);
                model.addConstr(cut1, GRB.GREATER_EQUAL,  - p_dt_double[d][t] + 0.5, "Gen1(" + (genCutCount) + ")");
                model.addConstr(cut1, GRB.LESS_EQUAL, params.M - p_dt_double[d][t] + 0.5, "Gen2(" + (genCutCount) + ")");
            }

        } catch (GRBException e) {
            e.printStackTrace();
        }
    }

    public SolStatus solve() {
        try {
            model.optimize();
            if (model.get(GRB.IntAttr.Status) == GRB.INFEASIBLE) {
                log.warn("Infeasible benders MP!");
                return SolStatus.INFEASIBLE;
            }

            objVal = model.get(GRB.DoubleAttr.ObjVal);
            queueLengh = model.get(GRB.DoubleAttr.X, L_dt);
            twdVal = model.get(GRB.DoubleAttr.X, twd);
            p = model.get(GRB.DoubleAttr.X, p_dt);
            cost_d_Val = model.get(GRB.DoubleAttr.X, cost_d);
            costOfAddedPhysicianValPerDay=model.get(GRB.DoubleAttr.X,costOfAddedPhysicianPerDay);

            for(int d=0;d<params.D;d++)
                 log.info("test2: {}", costOfAddedPhysicianValPerDay[d]);

            for(int d=0;d<7;d++)
                log.info("cost of d{} is {}",d,cost_d_Val[d]);

        } catch (GRBException e) {
            e.printStackTrace();
        }
        return SolStatus.OPTIMAL;
    }

    public double[][] getSol() throws GRBException {
        return p;
    }


    public void modifyTrustRegion(int iteration) throws GRBException {
        if (iteration > 1)
            removeAllTrustConstr();
        addAllTrustConstr();
    }


    private void removeAllTrustConstr() {
        for (int d = 0; d < params.D; d++) {
            try {
                model.remove(model.getConstrByName(String.format("trustRegion(%d)", d + 1)));
                for (int t = 0; t < params.T; t++) {
                    model.remove(model.getConstrByName(String.format("diff(%d,%d)", d + 1, t + 1)));
                    model.remove(model.getConstrByName(String.format("absDiffLeft(%d,%d)", d + 1, t + 1)));
                    model.remove(model.getConstrByName(String.format("absDiffRight(%d,%d)", d + 1, t + 1)));
                }
                for (int t = 0; t < params.T; t++) {
                    model.remove(model.getVarByName(String.format("absDiff[%d, %d]", d + 1, t + 1)));
                    model.remove(model.getVarByName(String.format("diff[%d,%d]", d + 1, t + 1)));
                }
            } catch (NullPointerException e) {
                log.warn("trust region constrs not found.");
            } catch (GRBException e) {
                e.printStackTrace();
            }
        }
    }

    // add trust region constrs
    private void addAllTrustConstr() throws GRBException {
        for (int d = 0; d < params.D; d++) {
            //绑定绝对值
            GRBVar[] absDiff = new GRBVar[params.T];
            GRBVar[] diff = new GRBVar[params.T];
            for (int t = 0; t < params.T; t++) {
                absDiff[t] = model.addVar(0, params.TRUST_SCOPE, 0, GRB.INTEGER, String.format("absDiff[%d, %d]", d + 1, t + 1));
                diff[t] = model.addVar(-params.TRUST_SCOPE, params.TRUST_SCOPE, 0, GRB.INTEGER, String.format("diff[%d,%d]", d + 1, t + 1));
            }
            for (int t = 0; t < params.T; t++) {
                GRBLinExpr expr = new GRBLinExpr();
                expr.addTerm(1, p_dt[d][t]);
                expr.addTerm(1, diff[t]);
                model.addConstr(expr, GRB.EQUAL, Math.round(p[d][t]), String.format("diff(%d,%d)", d + 1, t + 1));//上一轮的结果与本轮结果的比较

                model.addConstr(diff[t], GRB.LESS_EQUAL, absDiff[t], String.format("absDiffLeft(%d,%d)", d + 1, t + 1));
                GRBLinExpr minusDiff = new GRBLinExpr();
                minusDiff.addTerm(-1, diff[t]);
                model.addConstr(minusDiff, GRB.LESS_EQUAL, absDiff[t], String.format("absDiffRight(%d,%d)", d + 1, t + 1));
            }

            GRBLinExpr absSum = new GRBLinExpr();
            for (int t = 0; t < params.T; t++) {
                absSum.addTerm(1, absDiff[t]);
            }
            model.addConstr(absSum, GRB.LESS_EQUAL, params.TRUST_SCOPE, String.format("trustRegion(%d)", d + 1));
        }
    }

    public double getCost_d(int d){
        return cost_d_Val[d];
    }

}
