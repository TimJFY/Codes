package elementary_Data_Structures_Arrays_And_Strings;

import java.util.ArrayList;
import java.util.HashMap;

/* 
 * You are given a string, S, and a list of words, L, that are all of the same length. 
 * Find all starting indices of substring(s) in S that is a concatenation of each word in L exactly once and without any intervening characters.
 * 
 * For example, given:
 * S: "barfoothefoobarman"
 * L: ["foo", "bar"]
 * You should return the indices: [0,9]. (order does not matter).
 */

public class SubstrOfConcatenations {

	public SubstrOfConcatenations() {

	}

	private void addWord(String w, HashMap<String, Integer> words) {
		if (words.containsKey(w)) {
			words.put(w, words.get(w) + 1);
		} else {
			words.put(w, 1);
		}
	}

	private void removeWord(String w, HashMap<String, Integer> words) {
		if (!words.containsKey(w))
			return;
		if (words.get(w) > 1) {
			words.put(w, words.get(w) - 1);
		} else {
			words.remove(w);
		}
	}

	private int slideWindow(String S, int begin, int wordLen,
			HashMap<String, Integer> words) {
		// as the begin forwards, the valid word bounded by (begin, begin +
		// wordLen)
		// should be expected to find in the next check
		String old = S.substring(begin, begin + wordLen);
		addWord(old, words);
		return begin + wordLen;
	}

	public ArrayList<Integer> findSubstring(String S, String[] L) {
		ArrayList<Integer> indices = new ArrayList<Integer>();
		if (L.length == 0)
			return indices;

		int total = L.length, wordLen = L[0].length();

		// store the words and frequencies in L into a hash table
		HashMap<String, Integer> expectWords = new HashMap<String, Integer>();
		for (String w : L) {
			addWord(w, expectWords);
		}

		// find all concatenations
		// regard each word (length) as an interval, the possible begin should
		// in [0, length - 1]. 
		// from the overall perspective, all checks from length have been done when checking from 0
		// similarly, all checks from (length + 1) have been done when checking from 1... etc.
		for (int i = 0; i < wordLen; ++i) {
			// check if there are any concatenations
			// number of the current already matched valid words
			int count = 0;
			// all current expected words
			HashMap<String, Integer> collectWords = new HashMap<String, Integer>(
					expectWords);
			// the boundary should always leave sufficient space/length to match
			// the remaining words;
			// 'j' is the current index to check
			// 'begin' pivots the current valid substring
			for (int j = i, begin = i; j <= S.length() - (total - count)
					* wordLen
					&& begin <= S.length() - total * wordLen;) {
				String sub = S.substring(j, j + wordLen);
				if (!expectWords.containsKey(sub)) {
					// if encounter an invalid word, discard all the previous
					// matching and reset
					begin = j + wordLen;
					j = begin;
					count = 0;
					collectWords.putAll(expectWords);
				} else if (!collectWords.containsKey(sub)) {
					// if duplicate, keep forwarding begin by 1 word length
					begin = slideWindow(S, begin, wordLen, collectWords);
				} else {
					removeWord(sub, collectWords);
					j += wordLen;
					++count;
					// find a concatenation
					if (collectWords.isEmpty()) {
						indices.add(begin);
						begin = slideWindow(S, begin, wordLen, collectWords);
						--count;
					}
				}
			}
		}

		return indices;
	}

	public static void main(String[] args) {
		SubstrOfConcatenations substrOfConcatenations = new SubstrOfConcatenations();
		String s = "barfoothefoobarmanbarfodfoobafoobar";
		String[] l = { "foo", "bar" };
		ArrayList<Integer> results = substrOfConcatenations.findSubstring(s, l);
		System.out.print(results);
	}

}
