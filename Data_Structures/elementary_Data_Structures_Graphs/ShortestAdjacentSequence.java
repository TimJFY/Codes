package elementary_Data_Structures_Graphs;

import java.util.*;

public class ShortestAdjacentSequence {

	public int solution(int[] A) {
		if (A.length <= 2) {
			return A.length;
		}

		HashSet<Integer> visitedNodes = new HashSet<Integer>();
		Queue<Integer> queue = new LinkedList<Integer>();
		int minDepth = 2;
		boolean found = false;

		HashMap<Integer, HashSet<Integer>> neighboors = new HashMap<Integer, HashSet<Integer>>();
		for (int i = 0; i < A.length; i++) {
			if (neighboors.containsKey(A[i])) {
				if (i < A.length - 1) {
					neighboors.get(A[i]).add(A[i + 1]);
				}
				if (i > 0) {
					neighboors.get(A[i]).add(A[i - 1]);
				}
			} else {
				HashSet<Integer> nb = new HashSet<Integer>();
				if (i < A.length - 1) {
					nb.add(A[i + 1]);
				}
				if (i > 0) {
					nb.add(A[i - 1]);
				}
				neighboors.put(A[i], nb);
			}
		}

		queue.offer(A[0]);
		visitedNodes.add(A[0]);
		int preDepthCount = 1;
		while (!queue.isEmpty() && !found) {
			int curDepthCount = 0;
			while (preDepthCount > 0 && !found) {
				int curNode = queue.poll();
				HashSet<Integer> curNeighboors = neighboors.get(curNode);
				if (curNeighboors != null) {
					for (Integer curNB : curNeighboors) {
						if (curNB == A[A.length - 1]) {
							found = true;
							break;
						}
						if (!visitedNodes.contains(curNB)) {
							queue.offer(curNB);
							visitedNodes.add(curNB);
							curDepthCount++;
						}
					}
				}
				preDepthCount--;
			}
			preDepthCount = curDepthCount;
			if (!found) {
				minDepth++;
			}
		}
		return minDepth;
	}

	public static void main(String args[]) {
		ShortestAdjacentSequence sas = new ShortestAdjacentSequence();
		int A[] = new int[] { 2147483647, 2147483647, 2147483647 };
		System.out.println(sas.solution(A));
	}
}
