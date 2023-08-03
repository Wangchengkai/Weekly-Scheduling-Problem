package com.ben.benders_sp.branchNprice;

import com.ben.benders_sp.branchNprice.util.Coder;
import com.ben.benders_sp.branchNprice.util.Copier;
import com.ben.enums.SolStatus;
import com.ben.params;
import gurobi.*;
import lombok.Data;
import lombok.extern.slf4j.Slf4j;

import java.util.*;
import java.util.stream.Collectors;


import static com.ben.params.*;

/**
 * the master problem of column generation
 */

@Data
@Slf4j
public class CG_RMP {
	private Map<String, Integer> allScheduleSet;    //(Set W) The set of all schedules, where 'value' corresponds to the index in the 'schedules'.
	private List<List<Integer>> schedules;          //(Parameter a)The saved schedule
	private List<Integer> scheduleLengths;           //The on-duty time of each schedule

	//The set of optional schedules corresponding to (d,i),
	// where the first two dimensions are fixed and equal to d and i,
	// and the third dimension contains the indexes of the optional schedules in 'schedules'
	private List<Integer>[][] scheduleMap;
	private double objVal;                        //The objective value of RMP
	private double[][] p;                         // the expected number of physicians for each (d,i) from the master problem
	private boolean isIntSol = true;              // is the solution of RMP an integer solution?
	private double[] twd;                         // total working days of each physician
	private double costOfAddedPhysician;
	private double[] costOfAddedPhysicianPerDay;

	// gurobi related fields
	private GRBEnv env;
	private GRBModel model;
	private List<GRBVar>[][] y;

	private GRBConstr[][] minPhysicianConstrs;     // the minimum number of physicians (of each period) constraint
	private GRBConstr[][] oneRosterForEachConstrs; // each physician must choose one roster

	private double[][] minPhysicianDuals;          // the dual variables of the minimum number of physicians constraint
	private double[][] oneRosterForEachDuals;      //  the dual variables of the oneRosterForEach constraint

	private GRBVar[][] h;//医生上班

	/**
	 * Clone this model
	 *
	 * @return new model
	 */
	public CG_RMP clone() {
		return new CG_RMP(this);
	}

	/**
	 * initial this model from another RMP
	 *
	 * @param other another RMP
	 */
	public CG_RMP(CG_RMP other) {
		this.env = other.env;
		this.schedules = new ArrayList<>(other.schedules);
		this.scheduleLengths = new ArrayList<>(other.scheduleLengths);
		this.allScheduleSet = new HashMap<>(other.allScheduleSet);
		this.scheduleMap = Copier.copyArrayOfList(other.scheduleMap);
		this.p = other.p;
		this.minPhysicianDuals = new double[D][T];
		this.oneRosterForEachDuals = new double[N][D];

		//this.h=new
		try {
			this.model = new GRBModel(this.env);
			this.model.set(GRB.IntParam.OutputFlag, 0);
		} catch (GRBException e) {
			e.printStackTrace();
		}
	}

	public CG_RMP(GRBEnv env) {
		this.env = env;
		allScheduleSet = new HashMap<>();
		scheduleMap = new List[D][N];
		schedules = new ArrayList<>();
		scheduleLengths = new ArrayList<>();
		p = new double[D][T];
		minPhysicianDuals = new double[D][T];
		oneRosterForEachDuals = new double[N][D];
		try {
			model = new GRBModel(env);
			model.set(GRB.IntParam.OutputFlag, 0);
		} catch (GRBException e) {
			e.printStackTrace();
		}
	}

	// is integer solution
	public boolean isIntegerSol() {
		return isIntSol;
	}

	/**
	 *  update the minimum number of physicians constraint
	 *
	 * @param p_t the new array corresponding the minimum number constraint of physicians in each period
	 */
	public void resetP(double[][] p_t) {
		p = p_t;
		try {
			model.update();
			for (int d = 0; d < D; d++) {
				for (int t = 0; t < T; t++) {
					minPhysicianConstrs[d][t].set(GRB.DoubleAttr.RHS, p[d][t]);
				}
			}
		} catch (GRBException e) {
			e.printStackTrace();
		}
	}

	/**
	 * generate three initial schedules: empty schedule, night schedule, and a virtual full schedule
	 */
	public void addInitShifts() {

		int[][] initShifts = new int[3][T];
		for (int t = T - LEN_NIGHT_SCH; t < T; t++)
			initShifts[1][t] = 1;
		Arrays.fill(initShifts[2], 1);
		for (int i = 0; i < 3; i++) {
			List<Integer> sch = Arrays.stream(initShifts[i]).boxed().collect(Collectors.toList());
			allScheduleSet.put(Coder.encodeSchedule(sch), i);
			schedules.add(new ArrayList<>(sch));
			if (i == 0)
				scheduleLengths.add(0);
			else if (i == 1)
				scheduleLengths.add(LEN_NIGHT_SCH);
			else
				scheduleLengths.add(9999); // a large number
		}

		//  each physician can choose these three schedules
		List<Integer> indices = new ArrayList<>();
		for (int i = 0; i < schedules.size(); i++) indices.add(i);
		for (int i = 0; i < N; i++) {
			for (int d = 0; d < D; d++) {
				scheduleMap[d][i] = new ArrayList<>(indices);
			}
		}
	}

	//  add a column to the pool of each (d,i) pair
	public void addCol(int d, int i, List<Integer> schedule, int scheduleLength) throws GRBException {
		String encodedSchedule = Coder.encodeSchedule(schedule);
		if (!allScheduleSet.containsKey(encodedSchedule)) {
			allScheduleSet.put(encodedSchedule, schedules.size());
			schedules.add(new ArrayList<>(schedule));
			scheduleLengths.add(scheduleLength);
		}

		scheduleMap[d][i].add(allScheduleSet.get(encodedSchedule));
		double[] coef = new double[T];   // the coefficient of each schedule, corresponding to the 'a' in the model
		for (int t = 0; t < T; t++)
			coef[t] = 1.0 * schedule.get(t);

		model.update();
		GRBColumn col = new GRBColumn();
		col.addTerms(coef, minPhysicianConstrs[d]);
		col.addTerm(1.0, oneRosterForEachConstrs[i][d]);
		y[d][i].add(model.addVar(0, 1, scheduleLength, GRB.CONTINUOUS, col, String.format("y[%d,%d,%d]", d + 1, i + 1, scheduleMap[d][i].size())));
	}


	public SolStatus solve() throws GRBException {
		model.optimize();
		if (model.get(GRB.IntAttr.Status) == GRB.INFEASIBLE) {
			return SolStatus.INFEASIBLE;
		}

		objVal = model.get(GRB.DoubleAttr.ObjVal);

		// ensure that the solution is integer
		isIntSol = true;
		for (int d = 0; d < D; d++) {
			for (int i = 0; i < N; i++) {
				for (GRBVar var : y[d][i]) {
					double v = var.get(GRB.DoubleAttr.X);
					double vRounded = Math.round(v);
					if (Math.abs(v - vRounded) >= EPS) {
						isIntSol = false;
						break;
					}
				}
				if (!isIntSol) break;
			}
			if (!isIntSol) break;
		}

		for (int d = 0; d < D; d++) {
			for (int t = 0; t < T; t++)
				minPhysicianDuals[d][t] = minPhysicianConstrs[d][t].get(GRB.DoubleAttr.Pi);
		}

		for (int i = 0; i < N; i++) {
			for (int d = 0; d < D; d++) {
				oneRosterForEachDuals[i][d] = oneRosterForEachConstrs[i][d].get(GRB.DoubleAttr.Pi);
			}
		}

		// compute the total working hours of each day
		twd = new double[D];
		for (int d = 0; d < D; d++) {
			for (int t = 0; t < T; t++) {
				for (int i = 0; i < N; i++) {
					List<Integer> scheduleIndices = scheduleMap[d][i];
					for (int j = 0; j < scheduleIndices.size(); j++) {
						if (schedules.get(scheduleIndices.get(j)).get(t) == 1)
							twd[d] += y[d][i].get(j).get(GRB.DoubleAttr.X);
					}
				}
			}
		}

		//  compute the penalty cost of physician
		double[] test2=new double[IndexOfFinicPhysician];
		for (int m = 0; m < IndexOfFinicPhysician; m++) {
			List<Integer> scheduleIndices = scheduleMap[0][m];
			for(int i=0;i<scheduleIndices.size();i++){
				test2[m]+=y[0][m].get(i).get(GRB.DoubleAttr.X);
			}
			//log.info("test7 {} {} {}", 0, m, test2[m]);
		}


		costOfAddedPhysician=0;
		costOfAddedPhysicianPerDay=new double[D];
		for(int d=0;d<D;d++) {
			costOfAddedPhysicianPerDay[d]=0;
			for (int m = IndexOfFinicPhysician; m < N; m++) {
				costOfAddedPhysician += PuhishCostForAddedPhysician * (h[m][d].get(GRB.DoubleAttr.X));
				costOfAddedPhysicianPerDay[d] += PuhishCostForAddedPhysician* (h[m][d].get(GRB.DoubleAttr.X));
			}
		}

		return SolStatus.OPTIMAL;
	}

	public double[] getMinPhysicianDuals(int d) {
		return minPhysicianDuals[d];
	}

	public double getOneRosterForEachDual(int d, int i) {
		return oneRosterForEachDuals[i][d];
	}

	// return the branching node {d, i, t}
	public int[] findBranchPoint() throws GRBException {
		int[] branchingPoint = new int[3];
		double minDistToHalf = 1;  // q_di与0.5的最小距离
		for (int t = 0; t < T; t++) {
			boolean foundInT = false;
			minDistToHalf = 1;
			for (int i = 0; i < N; i++) {
				for (int d = 0; d < D; d++) {
					double q_di = 0;
					for (int j = 0; j < y[d][i].size(); j++) {
						if (schedules.get(scheduleMap[d][i].get(j)).get(t) == 1)
							q_di += y[d][i].get(j).get(GRB.DoubleAttr.X);
					}
					if (q_di > EPS && q_di < 1 - EPS) {
						foundInT = true;
						if (Math.abs(q_di - 0.5) < minDistToHalf - EPS) {
							minDistToHalf = Math.abs(q_di - 0.5);
							branchingPoint = new int[]{d, i, t};
						}
					}
				}
			}
			if (foundInT) break;
		}
		log.debug("branching point: (d,i,t)={},{},{}", branchingPoint[0], branchingPoint[1], branchingPoint[2]);
		log.debug("min distance to half: {}", minDistToHalf);
		return branchingPoint;
	}

	public double[][][] getSol() throws GRBException {
		double[][][] sol = new double[N][D][T];
		for (int i = 0; i < N; i++) {
			for (int d = 0; d < D; d++) {
				int roster = -1;
				for (int j = 0; j < y[d][i].size(); j++) {
					if (Math.abs(y[d][i].get(j).get(GRB.DoubleAttr.X) - 1) < EPS) {
						roster = scheduleMap[d][i].get(j);
						break;
					}
				}
				if(i>= IndexOfFinicPhysician &&roster<0)continue;;
				sol[i][d] = schedules.get(roster).stream().mapToDouble(k -> k).toArray();
			}
		}
		return sol;
	}

	public void initModel() throws GRBException {
		model.update();
		y = new ArrayList[D][N]; // y[d][i][j]: d for day, i for physician, j for shift
		minPhysicianConstrs = new GRBConstr[D][T];
		oneRosterForEachConstrs = new GRBConstr[N][D];


		// add variables
		for (int d = 0; d < D; d++) {
			for (int i = 0; i < N; i++) {
				y[d][i] = new ArrayList<>();
				List<Integer> scheduleIndices = scheduleMap[d][i];
				for (int j = 0; j < scheduleIndices.size(); j++) {
					// y relaxed to be continuous
					y[d][i].add(model.addVar(0, 1, scheduleLengths.get(scheduleIndices.get(j)), GRB.CONTINUOUS, String.format("y[%d,%d,%d]", d + 1, i + 1, j + 1)));
				}
			}
		}

		// the number of physicians who are on duty
		h=new GRBVar[N][D];
		for (int d = 0; d < D; d++) {
			for (int i = 0; i < N; i++) {
				h[i][d]= model.addVar(0, 1, PuhishCostForAddedPhysician, GRB.CONTINUOUS, String.format("h[%d,%d]", i+1,d + 1));
			}
		}

		// the length of y[d][i] === the length of scheduleMap[d][i]

		GRBVar z = model.addVar(0, 1, 0, GRB.CONTINUOUS, "z");

		// add constrs
		model.set(GRB.IntAttr.ModelSense, GRB.MINIMIZE);

		//the number of physicians in each period is depended on the schedule they choose.
		for (int d = 0; d < D; d++) {
			for (int t = 0; t < T; t++) {
				GRBLinExpr minPhysicianNum = new GRBLinExpr();
				for (int i = 0; i < N; i++) {
					List<Integer> scheduleIndices = scheduleMap[d][i];
					for (int j = 0; j < scheduleIndices.size(); j++) {
						if (schedules.get(scheduleIndices.get(j)).get(t) == 1)
							minPhysicianNum.addTerm(1.0, y[d][i].get(j));
					}
				}
				minPhysicianConstrs[d][t] = model.addConstr(minPhysicianNum, GRB.GREATER_EQUAL, p[d][t], String.format("minPhysicianNum(%d,%d)", d + 1, t + 1));

			}
		}

		// each physician must choose one and only one roster each day (changed)
		for (int d = 0; d < D; d++) {
			for (int i = 0; i < IndexOfFinicPhysician; i++) {
				GRBLinExpr rosters = new GRBLinExpr();
				List<Integer> scheduleIndices = scheduleMap[d][i];
				for (int j = 0; j < scheduleIndices.size(); j++)
					rosters.addTerm(1, y[d][i].get(j));

				oneRosterForEachConstrs[i][d] = model.addConstr(rosters, GRB.EQUAL, 1.0, String.format("OneRosterForEach(%d,%d)", i + 1, d + 1));
			}
		}

		// each physician must have at least one day off each week (changed)
		for (int i = 0; i < IndexOfFinicPhysician; i++) {
			GRBLinExpr oneDayOff = new GRBLinExpr();
			for (int d = 0; d < D; d++) {
				// 没有空排班可选时，就不加入，但如果所有7天都没有空排班可选，模型无解，通过一个额外变量z来控制。
				if (y[d][i].size() == 0 || scheduleMap[d][i].get(0) != 0) continue;
				oneDayOff.addTerm(1.0, y[d][i].get(0));
			}
			if (oneDayOff.size() == 0)
				model.addConstr(z, GRB.GREATER_EQUAL, 2, "infeasible");
			else
				model.addConstr(oneDayOff, GRB.GREATER_EQUAL, 1, String.format("OneDayOff(%d)", i + 1));
		}

		// have one day off after night shift
		for (int i = 0; i < IndexOfFinicPhysician; i++) {
			for (int d = 0; d < D - 1; d++) {
				if (!isNightShiftAvailable(scheduleMap[d][i]))
					continue;
				if (isEmptyShiftAvailable(scheduleMap[d + 1][i])) {
					if (scheduleMap[d][i].get(0) == 1)
						model.addConstr(y[d + 1][i].get(0), GRB.GREATER_EQUAL, y[d][i].get(0), String.format("NightOff(%d,%d)", d + 1, i + 1));
					else
						model.addConstr(y[d + 1][i].get(0), GRB.GREATER_EQUAL, y[d][i].get(1), String.format("NightOff(%d,%d)", d + 1, i + 1));
				} else {  // 如果下一天不可选休息，那么今天就不能选夜班
					if (scheduleMap[d][i].get(0) == 1)
						model.addConstr(y[d][i].get(0), GRB.LESS_EQUAL, 0, String.format("NightForbidden(%d,%d)", d + 1, i + 1));
					else
						model.addConstr(y[d][i].get(1), GRB.LESS_EQUAL, 0, String.format("NightForbidden(%d,%d)", d + 1, i + 1));
				}
			}
		}

		// the seconded physician (IndexOfFinicPhysician ~ N-1)
		// each seconded physician must choose one and only one roster each day (changed)
		for (int d = 0; d < D; d++) {
			for (int i = IndexOfFinicPhysician; i < N; i++) {
				GRBLinExpr rosters = new GRBLinExpr();
				List<Integer> scheduleIndices = scheduleMap[d][i];
				for (int j = 0; j < scheduleIndices.size(); j++)
					rosters.addTerm(1, y[d][i].get(j));
				rosters.addTerm(-1, h[i][d]);

				oneRosterForEachConstrs[i][d] = model.addConstr(rosters, GRB.EQUAL, 0.0, String.format("OneRosterForEach(%d,%d)", i + 1, d + 1));
			}
		}

		// the seconded physician is ordered by index
		for(int d=0;d<D;d++){
			for (int i = IndexOfFinicPhysician; i < N-1; i++) {
				GRBLinExpr physicianIndexDecreasing = new GRBLinExpr();
				physicianIndexDecreasing.addTerm(1.0, h[i][d]);
				physicianIndexDecreasing.addTerm(-1.0, h[i+1][d]);
				model.addConstr(physicianIndexDecreasing, GRB.GREATER_EQUAL, 0, String.format("h_decrease(%d)", i));
			}
		}

	}

	boolean isEmptyShiftAvailable(List<Integer> indices) {
		// if empty shift is available, it must be at index 0
		if (indices.isEmpty()) return false;
		return indices.get(0) == 0;
	}

	boolean isNightShiftAvailable(List<Integer> indices) {
		// if night shift is available, it must be at index 0 or 1
		if (indices.isEmpty()) return false;
		if (indices.get(0) == 1) return true;
		if (indices.size() > 1 && indices.get(1) == 1) return true;
		return false;
	}

	/**
	 * reset the scheduleMap according to the branching point
	 * @param branchingPoint branching point
	 * @param choice        0 means the branching point does not contain t
	 */
	public void resetScheduleMap(int[] branchingPoint, int choice) {
		int d = branchingPoint[0];
		int i = branchingPoint[1];
		int t = branchingPoint[2];

		List<Integer> indices = scheduleMap[d][i];
		if (choice == 0) { // the schedule containing t should be removed
			for (int j = indices.size() - 1; j >= 0; j--) {
				if (schedules.get(indices.get(j)).get(t) == 1)
					indices.remove(j);
			}
		} else { // the schedule not containing t should be removed
			for (int j = indices.size() - 1; j >= 0; j--) {
				if (schedules.get(indices.get(j)).get(t) != 1)
					indices.remove(j);
			}
		}
	}
}
