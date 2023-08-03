package com.ben.benders_sp;

import com.ben.benders_sp.branchNprice.BBNode;
import com.ben.benders_sp.branchNprice.util.Coder;
import com.ben.params;
import com.ben.enums.SolStatus;
import gurobi.GRBEnv;
import gurobi.GRBException;
import lombok.Data;
import lombok.extern.slf4j.Slf4j;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Stack;

import static com.ben.params.D;

/**
 * 分支定价框架
 */

@Slf4j
@Data
public class SPBranchNPrice {
	private double objVal;      //SP最优值
	private double[][][] incumbent;  // [i,d,t] 医生i在第d天的时段t是否上班
	private Stack<BBNode> nodes;
	private BBNode root;
	private int nodeCount;
	private int[][] staffing;  // 子问题求解得到的staffing
	private double[] twd;       // 每天的总工作时长
	private double costOfAddedPhysician;
	private double[] costOfAddedPhysicianPerDay;

	public SPBranchNPrice(GRBEnv env) throws GRBException {
		nodes = new Stack<>();
		nodeCount = 1;
		root = new BBNode(env, nodeCount);
	}

	/**
	 *	solve the column generation relaxation
	 * @return the status of the solution
	 * @throws GRBException
	 */
	public SolStatus solveLP(double[][] p_t) throws GRBException {
		resetStaffing(p_t);
		while (!nodes.empty()) {
			log.info("the number of nodes {}",nodes.size());
			BBNode node = nodes.pop();
			if (node.solve() == SolStatus.OPTIMAL) {  // status = 1 --> 最优解
				objVal = node.getObj();
				twd = node.getTWD();
				costOfAddedPhysician = node.getCostOfAddPhysicians();
				costOfAddedPhysicianPerDay=node.getCostOfAddphysiciansPerDay();
				//incumbent = node.getSol();//test
			} else {
				log.warn("Infeasible BP node!");
				return SolStatus.INFEASIBLE;
			}
		}
		//extractStaffing();//test
		return SolStatus.OPTIMAL;
	}

	/**
	 * Solve the branch and price
	 * @return  status of the solution
	 * @throws GRBException
	 */
	public SolStatus solve(double[][] p_t) throws GRBException {
		resetStaffing(p_t);
		boolean feasible = false;
		while (!nodes.empty()) {
			BBNode node = nodes.pop();
            log.info("============== solving node {} ==============", node.getID());
			if (node.solve() == SolStatus.OPTIMAL) {  // status = 1 --> 最优解
				feasible = true;
				if (node.getRmp().isIntegerSol()) {
					if (node.getObj() < objVal - params.EPS) {
						log.info("Found integer sol, current node obj val: {}", node.getObj());
						objVal = node.getObj();
						twd = node.getTWD();
						costOfAddedPhysician = node.getCostOfAddPhysicians();
						costOfAddedPhysicianPerDay=node.getCostOfAddphysiciansPerDay();
						incumbent = node.getSol();
						log.info("t1 {}", objVal);

					}
				} else if (node.getObj() < objVal -params.EPS) {
					log.info("t2 {}", node.getObj());
					branch(node);
				}



			} else {
				log.warn("Infeasible BP node!");
				continue;
			}
		}
		if (feasible) {
			extractStaffing();
			return SolStatus.OPTIMAL;
		}
		return SolStatus.INFEASIBLE;
	}

	private void extractStaffing() {
		staffing = new int[params.D][params.T];
		for (int d = 0; d < params.D; d++) {
			for (int t = 0; t < params.T; t++) {
				for (int i = 0; i < params.N; i++) {
					staffing[d][t] += (int) Math.round(incumbent[i][d][t]);
				}
			}
		}
	}

	public void saveResult() throws IOException {
		BufferedWriter br = new BufferedWriter(new FileWriter("data/result.txt", true));  // true for append
		br.write("\n*****************  SP RESULT  *****************\n");
		br.write("objective value: " + objVal + "\n");
		br.write("costOfPhysician: " + costOfAddedPhysician + "\n");
		br.write("number of physician in each period: ");
		for (int d = 0; d < params.D; d++) {
			for (int t = 0; t < params.T; t++) {
				if (staffing[d][t] == 0)
					log.error("sp staffing is zero! {d,t}={},{}", d, t);
				br.write(String.valueOf(staffing[d][t]) + ", ");
			}
			br.write("\n");
		}

		br.write("the schedule of each physician: \n");
		for (int i = 0; i < params.N; i++) {
			for (int d = 0; d < params.D; d++) {
				for (int t = 0; t < params.T; t++) {
					br.write(String.valueOf(incumbent[i][d][t]) + ", ");
				}
			}
			br.write("\n");
		}
		;
		br.write("the cost of seconded physician (in each day) : \n");
		for (int d = 0; d < params.D; d++) {
			br.write(String.valueOf(costOfAddedPhysicianPerDay[d]));
			br.write("\n");
		}

		br.write("\n\n");
		br.close();
	}

	private void resetStaffing(double[][] p_t) {
		objVal = Integer.MAX_VALUE;
		root.resetP(p_t);
		nodes.add(root);
	}

	// branch
	public void branch(BBNode node) throws GRBException {
		log.debug("branching at node: {}, get children {} and {}", node.getID(), nodeCount + 1, nodeCount + 2);
		int[] branchingPoint = node.findBranchPoint();
		BBNode leftNode = node.clone(nodeCount + 1);
		node = node.clone(nodeCount + 2);

		String code = Coder.encodeDIPair(branchingPoint[0], branchingPoint[1]);


		node.getMustExclude().get(code).add(branchingPoint[2]);
		node.getRmp().resetScheduleMap(branchingPoint, 0);
		node.getRmp().initModel();

		leftNode.getMustInclude().get(code).add(branchingPoint[2]);
		leftNode.getRmp().resetScheduleMap(branchingPoint, 1);
		leftNode.getRmp().initModel();

		nodes.add(node);
		nodes.add(leftNode);
		nodeCount += 2;
	}
}
