package Mathematics_Algorithms;

import java.util.ArrayList;

/*
 * Given an unsorted array of positive integers. Find the number of triangles that can be formed with three different array elements as three sides of triangles. 
 * For a triangle to be possible from 3 values, the sum of any two values (or sides) must be greater than the third value (or third side).
 * For example, if the input array is {4, 6, 3, 7}, the output should be 3. There are three triangles possible {3, 4, 6}, {4, 6, 7} and {3, 6, 7}. 
 * Note that {3, 4, 7} is not a possible triangle.
 * As another example, consider the array {10, 21, 22, 100, 101, 200, 300}. There can be 6 possible triangles: {10, 21, 22}, {21, 100, 101}, {22, 100, 101}, {10, 100, 101}, {100, 101, 200} and {101, 200, 300}
 */
public class FindPossibleTriangles {
	public ArrayList<ArrayList<Integer>> findTriangles(int[] edgeCadidates) {
		ArrayList<ArrayList<Integer>> results = new ArrayList<ArrayList<Integer>>();
		int count = 0;
		for (int i = 0; i < edgeCadidates.length - 2; i++) {
			int k = i + 2;
			for (int j = i + 1; j < edgeCadidates.length - 1; j++) {
				while (k < edgeCadidates.length
						&& edgeCadidates[i] + edgeCadidates[j] > edgeCadidates[k]) {
					if (j != k) {
						ArrayList<Integer> result = new ArrayList<Integer>();
						result.add(edgeCadidates[i]);
						result.add(edgeCadidates[j]);
						result.add(edgeCadidates[k]);
						results.add(result);
					}
					k++;
				}
				count += k - j - 1;
			}
		}
		System.out.println("Total Solutions:" + count);
		return results;
	}

	public static void main(String[] args) {
		FindPossibleTriangles findPossibleTriangles = new FindPossibleTriangles();
		ArrayList<ArrayList<Integer>> triangles = findPossibleTriangles
				.findTriangles(new int[] { 10, 21, 22, 100, 101, 200, 300 });
		System.out.println(triangles);
	}
}
