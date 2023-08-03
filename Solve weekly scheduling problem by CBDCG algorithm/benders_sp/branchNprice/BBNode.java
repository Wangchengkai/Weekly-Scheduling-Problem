package com.ben.benders_sp.branchNprice;

import com.ben.benders_sp.branchNprice.util.Copier;
import com.ben.params;
import com.ben.enums.SolStatus;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBVar;
import lombok.Data;
import lombok.extern.slf4j.Slf4j;
import com.ben.benders_sp.branchNprice.util.Coder;

import java.util.*;

/**
 * Nodes in the branch-and-bound tree
 */

@Data
@Slf4j
public class BBNode {
	private GRBEnv env;
	private CG_RMP rmp;
	private CG_SP sp;
	private int ID;
	private Map<String, Set<Integer>> mustInclude;  // (d,i)-->The time slots that must be included in the subproblem.
	private Map<String, Set<Integer>> mustExclude;  // (d,i)-->The time slots that must not be included in the subproblem.

	public BBNode clone(int id) throws GRBException {
		BBNode nodeCloned = null;
		try {
			nodeCloned = new BBNode(env, id);
		} catch (GRBException e) {
			e.printStackTrace();
		}
		nodeCloned.rmp = rmp.clone();
		nodeCloned.sp = sp.clone();
		nodeCloned.setMustInclude(Copier.copyMapOfSet(mustInclude));
		nodeCloned.setMustExclude(Copier.copyMapOfSet(mustExclude));
		return nodeCloned;
	}

	public BBNode(GRBEnv env, int id) throws GRBException {
		ID = id;
		this.env = env;
		mustInclude = new HashMap<>();
		mustExclude = new HashMap<>();
		for (int i = 0; i < params.N; i++) {
			for (int d = 0; d < params.D; d++) {
				mustInclude.put(Coder.encodeDIPair(d, i), new HashSet<>());
				mustExclude.put(Coder.encodeDIPair(d, i), new HashSet<>());
			}
		}
		rmp = new CG_RMP(env);
		// 初始化默认的几个排班
		rmp.addInitShifts();
		rmp.initModel();
		sp = new CG_SP();
	}

	public double[][][] getSol() throws GRBException {
		return rmp.getSol();
	}

	public double getObj() {
		return rmp.getObjVal();
	}

	public double[] getTWD() {
		return rmp.getTwd();
	}

	public double getCostOfAddPhysicians() {
		return rmp.getCostOfAddedPhysician();
	}

	public double[] getCostOfAddphysiciansPerDay(){return rmp.getCostOfAddedPhysicianPerDay();}

	public void resetP(double[][] p_t) {
		rmp.resetP(p_t);
	}

	public SolStatus solve() throws GRBException {
		int iterCount = 0;
		boolean canImprove;
		do {
			canImprove = false;
			iterCount++;
			log.debug("iteration: {} ...", iterCount);
			if ( rmp.solve()== SolStatus.INFEASIBLE) {
				log.warn("Infeasible RMP for column generation.");
				return SolStatus.INFEASIBLE;
			}
			log.debug("CG-RMP OBJ: {}", rmp.getObjVal());
			for (int i = 0; i < params.N; i++) {
				for (int d = 0; d < params.D; d++) {
					sp.clear();  // Clear the information of the previous solution iteration.
					sp.setPrices(rmp.getMinPhysicianDuals(d), rmp.getOneRosterForEachDual(d, i));
					String code = Coder.encodeDIPair(d, i);
					if (SolStatus.INFEASIBLE == sp.solve(mustInclude.get(code), mustExclude.get(code),i)) continue;
					log.debug("schedule: {}", sp.getSchedule().toString());
					log.debug("CG-SP ({},{}) obj: {}", d, i, sp.getObjVal() - rmp.getOneRosterForEachDual(d, i));
					if (sp.getObjVal() <= -params.EPS) {
						canImprove = true;//reduced cost is negative
						rmp.addCol(d, i, sp.getSchedule(), sp.getScheduleLen());
					}
				}
			}
		} while (canImprove);
		return SolStatus.OPTIMAL;
	}

	/**
	 * Find the branch node(d,i,t)
	 */
	public int[] findBranchPoint() throws GRBException {
		int[] branchPoint = rmp.findBranchPoint();
		return branchPoint;
	}
}
