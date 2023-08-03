package com.ben.benders_mp.util;

import com.ben.params;

public class Coder {

    public static String encodeCut(double[] p, int d, double lb) {
        StringBuilder sb = new StringBuilder();
        sb.append((int) Math.round(lb))
                .append(',')
                .append(d)
                .append(',');

        for (int t = 0; t < params.T; t++) {
            sb.append((int) Math.round(p[t]));
        }
        return sb.toString();
    }

    public static String encodeCut(double[] p, int d, double lb,double lb1) {
        StringBuilder sb = new StringBuilder();
        sb.append((int) Math.round(lb))
                .append(',')
                .append(d)
                .append(',')
                .append((int) Math.round(lb1))
                .append(',');

        for (int t = 0; t < params.T; t++) {
            sb.append((int) Math.round(p[t]));
        }
        return sb.toString();
    }

    public static String encodeCut(double[] p, double lb,double lb1, int d1) {
        StringBuilder sb = new StringBuilder();
        sb.append((int) Math.round(lb))
                .append(',')
                .append((int) Math.round(lb1))
                .append(',')
                .append(d1)
                .append(',');

        for (int t = 0; t < params.T; t++) {
            sb.append((int) Math.round(p[t]));
        }
        return sb.toString();
    }

    public static String encodeCutWeekly(double[][] p, double lb) {
        StringBuilder sb = new StringBuilder();
        sb.append((int) Math.round(lb))
                .append(',');
        for (int d = 0; d < params.D; d++){
            for (int t = 0; t < params.T; t++) {
                sb.append((int) Math.round(p[d][t]));
            }
        }
        return sb.toString();
    }

    public static String encodeStaffing(int day, double[] p) {
        StringBuilder sb = new StringBuilder();
        sb.append(day);
        sb.append(",");
        for (int t = 0; t < params.T - 1; t++) {
            sb.append((int) Math.round(p[t]));
            sb.append(",");
        }
        sb.append((int) Math.round(p[params.T - 1]));
        return sb.toString();
    }
}
