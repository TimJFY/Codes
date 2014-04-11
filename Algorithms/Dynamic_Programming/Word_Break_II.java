package Dynamic_Programming;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

public class Word_Break_II {
	public ArrayList<String> wordBreakII(String s, Set<String> dict) {
		int len = s.length();
		boolean f[] = new boolean[len + 1];
		ArrayList<ArrayList<Integer>> prev = new ArrayList<ArrayList<Integer>>();
		for (int i = 0; i <= len; i++) {
			prev.add(new ArrayList<Integer>());
		}
		f[0] = true;
		for (int i = 1; i <= len; i++) {
			for (int j = i - 1; j >= 0; j--) {
				if (f[j] && dict.contains(s.substring(j, i))) {
					prev.get(i).add(j);
					f[i] = true;
				}
			}
		}
		ArrayList<String> result = new ArrayList<String>();
		buildResult(s, prev, len, result, "");
		return result;
	}

	private void buildResult(String s, ArrayList<ArrayList<Integer>> prev,
			int end, ArrayList<String> result, String current) {
		ArrayList<Integer> prevs = prev.get(end);
		for (int i = 0; i < prevs.size(); i++) {
			int n = prevs.get(i);
			String sub = s.substring(n, end);
			String r = sub;
			if (!current.equals("")) {
				r = r + " " + current;
			}
			if (n == 0) {
				result.add(r);
			} else {
				buildResult(s, prev, n, result, r);
			}
		}
	}

	public static void main(String args[]) {

		Word_Break_II wb = new Word_Break_II();
		String s = "catsanddog";

		Set<String> dict = new HashSet<String>();
		dict.add("cat");
		dict.add("cats");
		dict.add("and");
		dict.add("sand");
		dict.add("dog");

		System.out.println(wb.wordBreakII(s, dict));

	}
}
