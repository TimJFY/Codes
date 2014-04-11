package Dynamic_Programming;

import java.util.Arrays;

public class Longest_Arithmetic_Progress {

	public int findLAS(int[] array) {
		int n = array.length;
		Arrays.sort(array);

		// matrix[i][j] means the longest arithmetic sequence which uses ith
		// number as the first element and jth number as the second element
		int[][] matrix = new int[n][n];
		int maxLength = 2;
		int maxFirst = 0;
		int maxSecond = 1;

		for (int i = 0; i < n; i++) {
			matrix[i][n - 1] = 2;
			// matrix[i][i] = 1;
		}

		for (int j = n - 2; j >= 1; j--) {
			int i = j - 1;
			int k = j + 1;
			while (i >= 0 && k <= n - 1) {
				if (array[i] + array[k] > 2 * array[j]) {
					matrix[i][j] = 2;
					i--;
				} else if (array[i] + array[k] < 2 * array[j]) {
					k++;
				}
				// if (array[i] + array[k] == 2 * array[j])
				else {
					matrix[i][j] = matrix[j][k] + 1;
					if (maxLength < matrix[i][j]) {
						maxLength = matrix[i][j];
						maxFirst = i;
						maxSecond = j;
					}
					i--;
					k++;
				}
			}

			while (i >= 0) {
				matrix[i][j] = 2;
				i--;
			}
		}
		System.out.println("1st element: " + array[maxFirst]
				+ "; common difference: " + (array[maxSecond] - array[maxFirst]));

		return maxLength;
	}

	public static void main(String[] args) {
		Longest_Arithmetic_Progress las = new Longest_Arithmetic_Progress();
		System.out.println("Total length = "
				+ las.findLAS(new int[] { 7, 3, 2, 9, 4, 5, 1, 0 }));
	}
}
