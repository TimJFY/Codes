package elementary_Data_Structures_Graphs;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.TreeMap;

// Can only generate partial paths
// The assumption here is that for each ladder inside the dictionary, 
// we only consider one of its pre-ladders and ignore others

public class WorldLadderII_Partial {
	public ArrayList<ArrayList<String>> findLadders(String start, String end,
			HashSet<String> dict) {
		Queue<String> actionQueue = new LinkedList<String>();
		Set<String> visitedWords = new HashSet<String>();
		Map<String, String> backtrackMap = new TreeMap<String, String>();
		int nextLayerSize = 0;
		int curLayerSize = 0;
		int wordsProcessed = 1;
		boolean stop = false;
		ArrayList<ArrayList<String>> results = new ArrayList<ArrayList<String>>();

		if (start == end) {
			if (dict.contains(start)) {
				ArrayList<String> oneKey = new ArrayList<String>();
				oneKey.add(start);
				oneKey.add(end);
				results.add(oneKey);
				return results;
			}
		}
		actionQueue.add(start);
		visitedWords.add(start);
		while (!actionQueue.isEmpty()) {
			if (wordsProcessed > curLayerSize) {
				if (stop) {
					break;
				}
				wordsProcessed = 1;
				curLayerSize = nextLayerSize;
				nextLayerSize = 0;
			}
			String w = actionQueue.poll();
			wordsProcessed++;
			for (String v : getLadders(w)) {
				if (v.equals(end)) {
					stop = true;
					ArrayList<String> oneKey = new ArrayList<String>();
					oneKey.add(v);
					while (w != null) {
						oneKey.add(w);
						w = backtrackMap.get(w);
					}
					results.add(reverse(oneKey));
				}

				else if (dict.contains(v)) {
					if (!visitedWords.contains(v)) {
						actionQueue.add(v);
						visitedWords.add(v);
						backtrackMap.put(v, w);
						nextLayerSize++;
					}
				}
			}

		}
		return results;
	}

	public Set<String> getLadders(String word) {
		Set<String> possibleLadders = new HashSet<String>();
		for (int i = 0; i < word.length(); i++) {
			char[] wordLetters = word.toCharArray();
			for (char c = 'a'; c <= 'z'; c++) {
				if (wordLetters[i] != c) {
					wordLetters[i] = c;
					possibleLadders.add(new String(wordLetters));
				}
			}
		}
		return possibleLadders;
	}

	public ArrayList<String> reverse(ArrayList<String> key) {
		ArrayList<String> reversedKey = new ArrayList<String>();
		for (int i = key.size() - 1; i >= 0; i--) {
			reversedKey.add(key.get(i));
		}
		return reversedKey;
	}

	public static void main(String[] args) {
		WorldLadderII_Partial worldLadderII = new WorldLadderII_Partial();
		String start = "red";
		String end = "tax";
		HashSet<String> dict = new HashSet<String>();
		dict.add("ted");
		dict.add("tex");
		dict.add("red");
		dict.add("tax");
		dict.add("tad");
		dict.add("den");
		dict.add("rex");
		dict.add("pee");

		System.out.println(worldLadderII.findLadders(start, end, dict));

	}
}
