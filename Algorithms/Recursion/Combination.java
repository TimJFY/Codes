package Recursion;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class Combination {

	static Set<List<Integer>> result = new HashSet<List<Integer>>();
	static ArrayList<Integer> set = new ArrayList<Integer>();
	static ArrayList<Boolean> used = new ArrayList<Boolean>();

	// with duplicates
	public void getAllCombinations2(ArrayList<Integer> curResult, int select,
			int start) {
		if (start == curResult.size() - 1) {
			// will generate permutations if use curResult directly
			List oneRslt = ((List) curResult.clone()).subList(0, select);

			result.add(oneRslt);
			return;
		}
		for (int i = start; i < curResult.size(); i++) {
			Collections.swap(curResult, start, i); 
			getAllCombinations2(curResult, select, start + 1);
			Collections.swap(curResult, i, start);
		}

	}

	public void getAllCombinations(ArrayList<Integer> curResult, int select,
			int start, int step) {
		if (step == select) {
			List oneRslt = curResult;
			result.add(oneRslt);
			return;
		}
		for (int i = start; i < set.size(); i++) {
			ArrayList<Integer> curR = new ArrayList<Integer>(curResult);
			curR.add(set.get(i));
			getAllCombinations(curR, select, i + 1, step + 1);
		}

	}

	public void getAllPermutations(ArrayList<Integer> curResult,
			int step) {
		if (step == set.size()) {
			List oneRslt = curResult;
			result.add(oneRslt);
			return;
		}

		for (int i = 0; i < set.size(); i++) {
			ArrayList<Integer> curR = new ArrayList<Integer>(curResult);
			if (!used.get(i)) {
				used.set(i, true);
				curR.add(set.get(i));
				getAllPermutations(curR, step + 1);
				//curR.remove(curR.size() - 1);
				used.set(i, false);
			}
		}

	}

	public static void main(String[] args) {
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		int n = 1;
		int k = 1;
		try {
			System.out.println("n:");
			n = Integer.parseInt(br.readLine());
			System.out.println("k:");
		    k = Integer.parseInt(br.readLine());
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		for (int i = 1; i <= n; i++) {
			set.add(i);
			used.add(false);
		}

		Combination c = new Combination();

		// ArrayList<Integer> pool = new ArrayList<Integer>();
		// for (int i = 1; i <= x; i++) {
		// pool.add(i);
		// }
		// c.getAllCombinations2(pool, y, 0);
		c.getAllCombinations(new ArrayList<Integer>(), k, 0, 0);
		System.out.println("C("+ n +" ,"+ k +") = " + result.size());
		System.out.println(result);
		System.out.println();
		result.clear();
		c.getAllPermutations(new ArrayList<Integer>(), 0);	
		System.out.println("A("+ n +" ,"+ n +") = " + result.size());
		System.out.println(result);
	}

}
