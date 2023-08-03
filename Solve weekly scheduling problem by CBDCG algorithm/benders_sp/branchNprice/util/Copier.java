package com.ben.benders_sp.branchNprice.util;

import java.util.*;

/**
 * @author: Ben
 * @create: 2021-07-18 12:56
 * @content:
 **/
public class Copier {

	public static <E, F> Map<E, List<F>> copyMapOfList(Map<E, List<F>> map) {
		Map<E, List<F>> copy = new HashMap<>();
		for (Map.Entry<E, List<F>> entry : map.entrySet()) {
			copy.put(entry.getKey(), new ArrayList<F>(entry.getValue()));
		}
		return copy;
	}

	public static <E, F> Map<E, Set<F>> copyMapOfSet(Map<E, Set<F>> map) {
		Map<E, Set<F>> copy = new HashMap<>();
		for (Map.Entry<E, Set<F>> entry : map.entrySet()) {
			copy.put(entry.getKey(), new HashSet<>(entry.getValue()));
		}
		return copy;
	}

	public static List<Integer>[][] copyArrayOfList(List<Integer>[][] arr) {
		List<Integer>[][] ret = new List[arr.length][arr[0].length];
		for (int d = 0; d < arr.length; d++) {
			for (int i = 0; i < arr[0].length; i++) {
				ret[d][i] = new ArrayList<>(arr[d][i]);
			}
		}
		return ret;
	}
}
