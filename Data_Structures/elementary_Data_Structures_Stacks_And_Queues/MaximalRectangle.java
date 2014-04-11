package elementary_Data_Structures_Stacks_And_Queues;

public class MaximalRectangle {

	public int compute(char[][] matrix) {

		int row = matrix.length;

		if (row == 0) {
			return 0;
		}

		int col = matrix[0].length;

		int[][] rowStat = new int[row][col];
		int[][] maxRecStat = new int[row][col];

		int result = 0;

		for (int i = 0; i < row; i++) {
			int consecutiveOnes = 1;
			for (int j = 0; j < col; j++) {
				if (matrix[i][j] == '1') {
					rowStat[i][j] = consecutiveOnes;
					consecutiveOnes++;
				} else {
					consecutiveOnes = 1;
				}
			}
		}

		showMatrix(rowStat);
		System.out.println();
		for (int k = 0; k < col; k++) {
			for (int l = 0; l < row; l++) {
				int minimal = Integer.MAX_VALUE;
				int validHight = 1;
				int start = l;
				while (start < row && rowStat[start][k] > 0) {
					minimal = Math.min(minimal, rowStat[start][k]);
					int curMax = validHight * minimal;
					if (maxRecStat[l][k] < curMax) {
						maxRecStat[l][k] = curMax;
						if (maxRecStat[l][k] > result) {
							result = maxRecStat[l][k];
						}
					}
					validHight++;
					start++;
				}
			}
		}

		showMatrix(maxRecStat);
		System.out.println();
		return result;
	}

	public void showMatrix(int[][] matrix) {
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				System.out.print(matrix[i][j] + " ");
			}
			System.out.println();
		}
	}

	public static void main(String[] args) {
		MaximalRectangle maximalRectangle = new MaximalRectangle();

		char[][] matrix = { { '0', '1', '1', '0', '1' },
				{ '1', '1', '0', '1', '0' }, { '0', '1', '1', '1', '0' },
				{ '1', '1', '1', '1', '0' }, { '1', '1', '1', '1', '1' },
				{ '0', '0', '0', '0', '0' } };

		System.out.println(maximalRectangle.compute(matrix));

	}
}
