package com.ben;

import com.ben.benders_mp.MP;
import com.ben.benders_sp.SPBranchNPrice;
import com.ben.enums.SolStatus;
import gurobi.GRBEnv;
import gurobi.GRBException;
import lombok.extern.slf4j.Slf4j;

import java.io.IOException;
import java.util.Arrays;

/**
 * Benders framework
 */

@Slf4j
public class Benders {
    private GRBEnv env;
    private MP mp;
    private SPBranchNPrice sp;
    private double logicBasedObj;   // logic based 最新的最优目标值

    public Benders(GRBEnv env) throws GRBException {
        this.env = env;
        mp = new MP(env);
    }

    /**
     * First phase : relaxed SP
     * Second phase: integer SP
     *
     * @throws GRBException
     * @throws IOException
     */
    public void twoPhaseLogicBased() throws GRBException, IOException {
        logicBasedObj = Double.MAX_VALUE;
        int iteration = 0;

        // stage 1 with linear SP
        log.info("iterate with linear bender sp");
        do {
            iteration++;
            log.info("\n -------- MP: iteration {} --------\n", iteration);
            if (mp.solve() == SolStatus.INFEASIBLE) break;
            log.info("mp obj: {}", mp.getObjVal());
            sp = new SPBranchNPrice(env);       // re-initial the SP
            mp.saveResult();

            log.info("\n -------- SP --------\n");
            if (SolStatus.INFEASIBLE == sp.solveLP(mp.getSol())) {
                log.error("no solution to benders SP!");
                double[][] sol = mp.getSol();
                for (double[] d : sol) {
                    System.out.println(Arrays.toString(d));
                }
                break;
            }
            log.info("sp obj: {}", sp.getObjVal());

            //due to the consideration of penalty cost:
            double currentObjOfSP = 0;
            for(int d=0;d<params.D;d++)
                currentObjOfSP += mp.getCost_d(d);// add penalty cost
            currentObjOfSP += sp.getObjVal();// add schedule time
            log.info("sp obj(include punish): {}", currentObjOfSP);
            log.info("sp obj(include punish): {}", sp.getCostOfAddedPhysician());

            if (currentObjOfSP < logicBasedObj) {
                // update the upper bound of this node
                logicBasedObj = currentObjOfSP;
            }
            if (mp.getObjVal() > currentObjOfSP - params.EPS) {
                System.out.println("optimal val: " + logicBasedObj);
                break;
            }

            //cut
            mp.addCut(sp.getTwd(),sp.getCostOfAddedPhysicianPerDay());
            log.info("test3:{}", sp.getCostOfAddedPhysicianPerDay()[1]);

            mp.modifyTrustRegion(iteration);

        } while (true);

        // stage 2 with integer SP
        log.info("iterate with integer bender sp");
        iteration = 0;
        logicBasedObj = Double.MAX_VALUE;
        do {
            iteration++;
            log.info("\n -------- MP: iteration {} --------\n", iteration);
            if (mp.solve() == SolStatus.INFEASIBLE) break;
            log.info("mp obj: {}", mp.getObjVal());
            sp = new SPBranchNPrice(env);       // 重新初始化SP
            log.info("\n -------- SP --------\n");

            if (SolStatus.INFEASIBLE == sp.solve(mp.getSol())) {
                log.error("no solution to benders SP!");
                break;
            }
            log.info("sp obj: {}", sp.getObjVal());

            //因为考虑了惩罚成本
            double currentObjOfSP = 0;
            for(int d=0;d<params.D;d++)
                currentObjOfSP += mp.getCost_d(d);
            currentObjOfSP += sp.getObjVal();
            log.info("sp obj(include punish): {}", currentObjOfSP);

            if (currentObjOfSP < logicBasedObj) {    // 更新本节点的上界
                logicBasedObj = currentObjOfSP;
            }
            saveResult(mp, sp);
            if (mp.getObjVal() > currentObjOfSP - params.EPS) {
            //if (Math.abs(mp.getObjVal() - currentObjOfSP) <= params.EPS) {
                System.out.println("optimal val: " + logicBasedObj);
                return;
            }
            mp.addCut(sp.getTwd(),sp.getCostOfAddedPhysicianPerDay());
            mp.modifyTrustRegion(iteration);
        } while (true);
    }

    private void saveResult(MP mp, SPBranchNPrice sp) throws IOException {
        log.info("Generated " + mp.getGenCutCount() + " cuts in total.");
        mp.saveResult();
        sp.saveResult();
    }

    //	程序运行入口
    public static void main(String[] args) throws GRBException, IOException {
        GRBEnv env = new GRBEnv();
        Benders benders = new Benders(env);
        benders.twoPhaseLogicBased();
    }
}
