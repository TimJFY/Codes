package elementary_Data_Structures_Arrays_And_Strings;

import java.util.ArrayList;
import java.util.HashMap;

public class LongestUniqueSubstring {

	// To store the start indices of all Longest Unique Substrings
	public ArrayList<Integer> startIndices;

	public int findLUSubstring(String s) {
		// To store the result
		int maxLength = 0;
		if (s.length() == 0 || s == null) {
			return maxLength;
		}

		startIndices = new ArrayList<Integer>();
		int curStartIndex = 0;
		// To store the length of current substring
		int curLength = 0;
		HashMap<Character, Integer> visited = new HashMap<Character, Integer>();

		for (int i = 0; i < s.length(); i++) {
			char c = s.charAt(i);
			/*
			 * If the current character never presents before or it is not part
			 * of the current NRCS(Non-Repeating Character Substring), then
			 * expand the window of NRCS directly
			 */
			if (!visited.containsKey(c) || i - visited.get(c) > curLength) {
				curLength++;
			}
			/*
			 * If the current character is present in currently considered NRCS,
			 * then update NRCS to start from the next character of previous
			 * instance.
			 */
			else {
				curStartIndex = visited.get(c) + 1;
				curLength = i - curStartIndex + 1;
			}
			// update the index of current character
			visited.put(c, i);

			if (curLength >= maxLength) {
				if (curLength > maxLength) {
					startIndices.clear();
				}
				maxLength = curLength;
				startIndices.add(curStartIndex);
			}
		}
		return maxLength;
	}

	public void showSubstrings(String s, int maxLen) {
		for (Integer i : startIndices) {
			System.out.println(s.substring(i, i + maxLen));
		}
	}

	public static void main(String[] args) {
		LongestUniqueSubstring longestUniqueSubstring = new LongestUniqueSubstring();
		String s = "GEEKSFORGEEKS";
		int maxLen = longestUniqueSubstring.findLUSubstring(s);
		System.out.println("maxLen = " + maxLen);
		System.out.println();
		longestUniqueSubstring.showSubstrings(s, maxLen);
	}
}
