package Recursion;

/**
 * LeetCode - Minimum Window Substring
 * 
 * Given a string S and a string T, find the minimum window in S which will
 * contain all the characters in T in complexity O(n).
 * 
 * For example, S = "ADOBECODEBANC" T = "ABC" Minimum window is "BANC".
 * 
 * Note: If there is no such window in S that covers all characters in T, return
 * the emtpy string "".
 * 
 * If there are multiple such windows, you are guaranteed that there will always
 * be only one unique minimum window in S.
 */

public class MinimumWindowSubSring {
	public String minWindow(String S, String T) {
		// Note: The Solution object is instantiated only once and is reused by
		// each test case.
		char[] chS = S.toCharArray(), chT = T.toCharArray();
		int[] sum = new int[256];
		boolean[] flag = new boolean[256];
		int count = 0;
		for (int i = 0; i < chT.length; i++) {
			sum[chT[i]]++;
			flag[chT[i]] = true;
			count += sum[chT[i]] == 1 ? 1 : 0;
		}
		String answer = "";
		// the left and right boundary of sliding window
		int left = 0, right = 0;
		while (true) {
			// slide the right boundary to find the first solution window -
			// complete match
			while (count > 0 && right < chS.length) {
				sum[chS[right]] -= flag[chS[right]] ? 1 : 0;
				count -= sum[chS[right]] == 0 && flag[chS[right]] ? 1 : 0;
				right++;
			}
			// reach the end of S, but still no solution is detected
			if (count > 0) {
				break;
			}
			// slide the left boundary until encounter the first 'overdraw'
			// character
			// Note: to pivot not only the first but also the 'overdraw'
			// character is to slide the boundary as left as possible to ensure
			// the optimal window size, because if the current character is a
			// 'hit' but not 'overdrew', you can always use the remaining
			// characters to accomplish a full match
			while (left < chS.length && count == 0) {
				if ("".equals(answer) || right - left < answer.length()) {
					answer = S.substring(left, right);
				}
				sum[chS[left]] += flag[chS[left]] ? 1 : 0;
				count += sum[chS[left]] == 1 && flag[chS[left]] ? 1 : 0;
				left++;
			}
		}
		return answer;
	}

	public static void main(String[] args) {
		MinimumWindowSubSring mws = new MinimumWindowSubSring();
		System.out.println(mws.minWindow("ADOBECODEBANC", "ABC"));
	}
}
