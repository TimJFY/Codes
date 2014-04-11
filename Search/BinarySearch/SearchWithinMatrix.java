package BinarySearch;

import java.util.ArrayList;

//Wrong!!!
public class SearchWithinMatrix {
	public int[] searchBasedOnDiagonal(int[][] matrix, int target) {
		int[] position = new int[] { -1, -1 };
		int lengthX = matrix.length;
		int lengthY = matrix[0].length;
		// divide into squares
		int longSide = Math.max(lengthX, lengthY);
		int shortSide = Math.min(lengthX, lengthY);
		int pieces = longSide / shortSide;
		int remainder = longSide - pieces * shortSide;

		// binary search in the squares

		for (int i = 0; i < pieces; i++) {
			int start = i * shortSide;
			int end = i * shortSide + shortSide - 1;
			while (start < end) {
				int mid = (start + end) / 2;
				if (matrix[mid][mid] == target) {
					position[0] = position[1] = mid;
					return position;
				}
				if (matrix[mid][mid] > target) {
					end = mid - 1;
				} else {
					start = mid + 1;
				}
			}

			if (start == end) {
				// check the 9 grids, 
				// in fact, 9 grids is not enough, so this approach is invalid
				ArrayList<int[]> grids = new ArrayList<int[]>();
				grids.add(new int[] { matrix[start][end], start, end });
				if (start > 0) {
					grids.add(new int[] { matrix[start - 1][end], start - 1,
							end });
					grids.add(new int[] { matrix[start - 1][end - 1],
							start - 1, end - 1 });
					grids.add(new int[] { matrix[start][end - 1], start,
							end - 1 });
					if (start < shortSide - 1) {
						grids.add(new int[] { matrix[start + 1][end],
								start + 1, end });
						grids.add(new int[] { matrix[start + 1][end + 1],
								start + 1, end + 1 });
						grids.add(new int[] { matrix[start - 1][end + 1],
								start - 1, end + 1 });
						grids.add(new int[] { matrix[start][end + 1], start,
								end + 1 });
						grids.add(new int[] { matrix[start + 1][end - 1],
								start + 1, end - 1 });

					}
				} else {
					grids.add(new int[] { matrix[start + 1][end], start + 1,
							end });
					grids.add(new int[] { matrix[start + 1][end + 1],
							start + 1, end + 1 });
					grids.add(new int[] { matrix[start][end + 1],
							start - 1, end + 1 });
				}

				for (int[] value : grids) {
					if (value[0] == target) {
						position[0] = value[1];
						position[1] = value[2];
						return position;
					}
				}
			}
		}

		// binary search in the remainders(by rows or by columns)

		return position;
	}

	public static void main(String[] args) {
		int[][] matrix = new int[][] { { 10, 20, 30 }, { 11, 21, 31 },
				{ 12, 22, 32 } };
		int result[];
		SearchWithinMatrix searchWithinMatrix = new SearchWithinMatrix();
		result = searchWithinMatrix.searchBasedOnDiagonal(matrix, 12);
		System.out.println("Location: X = " + result[0] + ", Y = " + result[1]);

	}
}
