package com.ben.benders_sp.branchNprice.util;

import java.util.List;

/**
 * @author: Ben
 * @create: 2021-07-11 19:14
 * @content:
 **/
public class Coder {
    public static String encodeDIPair(int d, int i) {
        StringBuilder sb = new StringBuilder();
        sb.append(d);
        sb.append("#");
        sb.append(i);
        return sb.toString();
    }

    public static String encodeVWN(int v, int w, int n) {
        StringBuilder sb = new StringBuilder();
        sb.append(v);
        sb.append("#");
        sb.append(w);
        sb.append("#");
        sb.append(n);
        return sb.toString();
    }

    public static String encodeVWNR(int v, int w, int n,int r) {
        StringBuilder sb = new StringBuilder();
        sb.append(v);
        sb.append("#");
        sb.append(w);
        sb.append("#");
        sb.append(n);
        sb.append("#");
        sb.append(r);
        return sb.toString();
    }

    public static String encodeSchedule(List<Integer> schedule) {
        StringBuilder sb = new StringBuilder();
        for (Integer i : schedule) sb.append(i);
        return sb.toString();
    }
}
