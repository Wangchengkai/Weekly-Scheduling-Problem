package com.ben.benders_sp.branchNprice;

import com.ben.benders_sp.branchNprice.util.Coder;
import com.ben.params;
import com.ben.enums.SolStatus;
import lombok.Data;
import lombok.extern.slf4j.Slf4j;

import java.util.*;

import static java.lang.Math.min;

/**
 * the sub problem of column generation
 */

@Data
@Slf4j
public class CG_SP {
    private double commonDualVals;         // common dual values from RMP
    private double[] dualVals;             // dual values from RMP
    private double[] modDualVals;          // the period that must be included is set to M, and the period that must be excluded is set to -M
    private double objVal;                 // the optimal value of SP
    private List<Integer> schedule;        // the generated column, which is also the optimal solution of SP
    private int scheduleLen;               // the time length of the generated column
    private int maxT;
    private HashMap<String, Double> objMemo;          // the memo of the optimal value of dynamic programming, key is coded (v,w,n), val is the minimum target value
    private HashMap<String, List<Integer>> solMemo;   // the memo of the optimal solution of dynamic programming, key is coded (v,w,n), val is the period included in the optimal solution
    private Set<Integer> mustInclude;
    private Set<Integer> mustExclude;


    public CG_SP() {
        maxT = params.T - 8;
        objMemo = new HashMap<>();
        solMemo = new HashMap<>();
    }

    public CG_SP clone() {
        return new CG_SP();
    }

    public void setPrices(double[] dualVals, double commonDualVals) {
        this.dualVals = dualVals;
        this.commonDualVals = commonDualVals;
    }

    /**
     * reflect the period that must be included and the period that must be excluded to modDualVals
     * @param mustInclude  the period that must be included
     * @param mustExclude  the period that must be excluded
     */
    private void setConstraints(final Set<Integer> mustInclude, final Set<Integer> mustExclude) {
        this.mustInclude = mustInclude;
        this.mustExclude = mustExclude;
        this.modDualVals = Arrays.copyOf(dualVals, dualVals.length);
        for (Integer t : mustInclude)
            modDualVals[t] = params.M;
        for (Integer t : mustExclude)
            modDualVals[t] = -params.M;
    }

    public SolStatus solve(Set<Integer> mustInclude, Set<Integer> mustExclude,int i) {
        HashSet<Integer> result = new HashSet<>(mustExclude);
        result.retainAll(mustInclude);
        if (result.size() > 0) return SolStatus.INFEASIBLE;

        setConstraints(mustInclude, mustExclude);
        solveWithDP(i);
        return checkSchedule();
    }

    // there may be a situation where the result cannot satisfy mustInclude and mustExclude, in this case, -1 should be returned
    private SolStatus checkSchedule() {
        for (Integer t : mustExclude)
            if (schedule.get(t) == 1) return SolStatus.INFEASIBLE;

        for (Integer t : mustInclude)
            if (schedule.get(t) == 0) return SolStatus.INFEASIBLE;
        return SolStatus.OPTIMAL;
    }

    public void solveWithDP(int i) {
        // ensure the minimum value
        int r1=-1;int r2=-1;
        double retOne = Double.MAX_VALUE;
        double retTwo = Double.MAX_VALUE;

        for(int r=1;r<= params.S_MaxForOnePhy;r++){
            if( dp(0, maxT - 1, 1, r)<retOne) {
                retOne=dp(0, maxT - 1, 1, r);
                r1=r;
            }
        }
        for(int r=1;r<= params.S_MaxForOnePhy;r++){
            if( dp(0, maxT - 1, 2, r)<retTwo) {
                retTwo=dp(0, maxT - 1, 2, r);
                r2=r;
            }
        }

        // get the true objective value and the corresponding schedule
        List<Integer> sol;
        if (retOne > retTwo) {
            sol = solMemo.get(Coder.encodeVWNR(0, maxT - 1, 2, r2));
        } else {
            sol = solMemo.get(Coder.encodeVWNR(0, maxT - 1, 1, r1));
        }

        if(i<params.IndexOfFinicPhysician)objVal = -commonDualVals;
        else objVal=-commonDualVals+params.PuhishCostForAddedPhysician;

        schedule = new ArrayList<>(params.T);
        for (int t = 0; t < params.T; t++) schedule.add(0);
        for (Integer t : sol) {
            schedule.set(t, 1);
            objVal += 1 - dualVals[t];
        }

        scheduleLen = sol.size();
    }

    /**
     * dynamic programming
     * @param v start period
     * @param w end period (in the forward recursion, w is a fix period)
     * @param n number of shifts
     * @param r total number of on-duty periods
    **/
    private double dp(int v, int w, int n, int r) {
        if (objMemo.containsKey(Coder.encodeVWNR(v, w, n, r))) return objMemo.get(Coder.encodeVWNR(v, w, n, r));
        if (n == 0 || w - v + 1 <= params.S_LB || r < params.S_LB) {
            objMemo.put(Coder.encodeVWNR(v, w, n,r), 0.0);
            solMemo.put(Coder.encodeVWNR(v, w, n,r), new ArrayList<>());
            return 0;
        }

        double ret = Double.MAX_VALUE;
        List<Integer> optSol = new ArrayList<>();

        // do not choose v
        if (dp(v + 1, w, n, r) < ret) {
            String key = Coder.encodeVWNR(v + 1, w, n,r);
            ret = objMemo.get(key);
            optSol = solMemo.get(key);
        }

        // choose v
        int optD = -1;
        for (int d = params.S_LB; d <= min(params.S_UB,r); d++) {
            double cur = dp(v + d + params.G, w, n - 1, r-d);
            for (int t = v; t < min(v + d,params.T); t++)
                cur += (1 - modDualVals[t]);
            if (cur < ret) {
                ret = cur;
                optD = d;
            }
        }

        if (optD != -1) {
            optSol = new ArrayList<>();
            for (int t = v; t < v + optD; t++)
                optSol.add(t);
            optSol.addAll(solMemo.getOrDefault(Coder.encodeVWNR(v + optD + params.G, w, n - 1, r-optD), new ArrayList<>()));
        }
        solMemo.put(Coder.encodeVWNR(v, w, n,r), new ArrayList<>(optSol));
        objMemo.put(Coder.encodeVWNR(v, w, n,r), ret);
        return ret;
    }

    public void clear() {
        objMemo = new HashMap<>();
        solMemo = new HashMap<>();
    }
}
