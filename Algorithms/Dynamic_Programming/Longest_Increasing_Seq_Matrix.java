package Dynamic_Programming;

/*
 *  find longest increasing sub sequence in 2d array. 
 *  (bit more expl..) 
 *  ex: finding length of the snake in snake game 
 *  --------- 
 *  the sequence must not be diagonally. 
 *  but it can be any like top-bootm,bottom-left-top ........ 
 *  increasing means one step 
 *  ex: 10,11,12,13 (correct) 12,14,15,20(wrong) 
 *  Ex: input: consider 4x4 grid 
 *  	2   3   4   5 
 *      4   5   10  11 
 *     	20  6   9   12 
 *      6   7   8   40
 *      output : 4 5 6 7 8 9 10 11 12
 * 
 */

public class Longest_Increasing_Seq_Matrix {
	private int originalMatrix[][];
	// create a new 2d array (let say named temp[][]) of same dimension as
	// original.
	// the element temp[i][j] contains the length of the LIS starting at
	// grid[i][j].
	private int tmepLISMatrix[][];

	private int x_max;
	private int y_max;

	public void searchLongestIncreasingSeq(int matrix[][]) {
		// start grid of the LIS
		int lis_start_x = 0;
		int lis_start_y = 0;
		// length of the LIS
		int length = 0;

		x_max = matrix.length;
		if (x_max == 0) {
			return;
		}
		y_max = matrix[0].length;

		tmepLISMatrix = new int[x_max][y_max];
		originalMatrix = matrix;

		// fill up the grid at tmepLISMatrix[i][j]
		// return the value of tmepLISMatrix[i][j], Dynamic Programming, bottom-up
		for (int i = 0; i < x_max; i++) {
			for (int j = 0; j < y_max; j++) {
				if (tmepLISMatrix[i][j] == 0) {
					this.fill(tmepLISMatrix, i, j);
				}
			}
		}

		// find the ending grid of the LIS
		for (int i = 0; i < x_max; i++) {
			for (int j = 0; j < y_max; j++) {
				if (tmepLISMatrix[i][j] > length) {
					lis_start_x = i;
					lis_start_y = j;
					length = tmepLISMatrix[i][j];
				}
			}
		}

		this.showLISSeq(lis_start_x, lis_start_y);
		this.showLISSeq2(lis_start_x, lis_start_y);
	}

	// for this problem, we only need to check 2 rules: boundary limit and logic condition,
	// because the monotonicity of the logic condition ensures that no grid can be visited more than once.
	// but if the logic condition is not monotonic, we need to use global tag variables to avoid repetitive visits  
	public int fill(int lisMatrix[][], int x, int y) {
		if (x < 0 || y < 0 || x >= x_max || y >= y_max)
			return 0;
		if (lisMatrix[x][y] != 0)
			return lisMatrix[x][y];

		int left = 0;
		int right = 0;
		int up = 0;
		int down = 0;

		// go left
		if (x > 0 && originalMatrix[x][y] == originalMatrix[x - 1][y] - 1) {
			left = fill(lisMatrix, x - 1, y);
		}
		// go right
		if (x < x_max - 1
				&& originalMatrix[x][y] == originalMatrix[x + 1][y] - 1) {
			right = fill(lisMatrix, x + 1, y);
		}
		// go down
		if (y > 0 && originalMatrix[x][y] == originalMatrix[x][y - 1] - 1) {
			right = fill(lisMatrix, x, y - 1);
		}
		// go up
		if (y < y_max - 1
				&& originalMatrix[x][y] == originalMatrix[x][y + 1] - 1) {
			right = fill(lisMatrix, x, y + 1);
		}

		int rowMax = Math.max(left, right);
		int columnMax = Math.max(up, down);
		lisMatrix[x][y] = Math.max(rowMax, columnMax) + 1;

		return lisMatrix[x][y];
	}

	public void showLISSeq(int x, int y) {
		int i = x;
		int j = y;

		System.out.print("Start from " + "[" + i + ", " + j + "]" + " = "
				+ originalMatrix[i][j] + "... ");

		while (tmepLISMatrix[i][j] != 1) {
			// up
			if (i > 0) {
				if (tmepLISMatrix[i - 1][j] == tmepLISMatrix[i][j] - 1
						&& originalMatrix[i - 1][j] == originalMatrix[i][j] + 1) {
					i = i - 1;
					System.out.print("[" + i + ", " + j + "]" + " = "
							+ originalMatrix[i][j] + "; ");
					continue;
				}
			}
			// down
			if (i < x_max - 1) {
				if (tmepLISMatrix[i + 1][j] == tmepLISMatrix[i][j] - 1
						&& originalMatrix[i + 1][j] == originalMatrix[i][j] + 1) {
					i = i + 1;
					System.out.print("[" + i + ", " + j + "]" + " = "
							+ originalMatrix[i][j] + "; ");
					continue;
				}
			}
			// left
			if (j > 0) {
				if (tmepLISMatrix[i][j - 1] == tmepLISMatrix[i][j] - 1
						&& originalMatrix[i][j - 1] == originalMatrix[i][j] + 1) {
					j = j - 1;
					System.out.print("[" + i + ", " + j + "]" + " = "
							+ originalMatrix[i][j] + "; ");
					continue;
				}
			}
			// right
			if (j < y_max - 1) {
				if (tmepLISMatrix[i][j + 1] == tmepLISMatrix[i][j] - 1
						&& originalMatrix[i][j + 1] == originalMatrix[i][j] + 1) {
					j = j + 1;
					System.out.print("[" + i + ", " + j + "]" + " = "
							+ originalMatrix[i][j] + "; ");
					continue;
				}
			}
		}
		System.out.println();
	}

	// if we need to print out all LISs, then use one ArrayList to store each path
	// and do clone() at each branch, similar to find all paths(between 2 nodes) in a graph 
	public void showLISSeq2(int x, int y) {
		int i = x;
		int j = y;
		System.out.print("[" + i + ", " + j + "]" + " = "
				+ originalMatrix[i][j] + "; ");
		if (tmepLISMatrix[i][j] != 1) {
			if (i > 0) {
				if (tmepLISMatrix[i - 1][j] == tmepLISMatrix[i][j] - 1
						&& originalMatrix[i - 1][j] == originalMatrix[i][j] + 1) {
					showLISSeq2(i - 1, j);

				}
			}
			if (i < x_max - 1) {
				if (tmepLISMatrix[i + 1][j] == tmepLISMatrix[i][j] - 1
						&& originalMatrix[i + 1][j] == originalMatrix[i][j] + 1) {
					showLISSeq2(i + 1, j);
				}
			}
			if (j > 0) {
				if (tmepLISMatrix[i][j - 1] == tmepLISMatrix[i][j] - 1
						&& originalMatrix[i][j - 1] == originalMatrix[i][j] + 1) {
					showLISSeq2(i, j - 1);
				}
			}

			if (j < y_max - 1) {
				if (tmepLISMatrix[i][j + 1] == tmepLISMatrix[i][j] - 1
						&& originalMatrix[i][j + 1] == originalMatrix[i][j] + 1) {
					showLISSeq2(i, j + 1);
				}
			}
		}
	}

	public static void main(String[] args) {
		int[][] origMatrix = new int[][] { { 1, 2, 3, 4 }, { 8, 3, 8, 5 },
				{ 8, 4, 5, 6 }, { 0, 0, 8, 7 } };
		Longest_Increasing_Seq_Matrix longest_Increasing_Seq_Matrix = new Longest_Increasing_Seq_Matrix();
		longest_Increasing_Seq_Matrix.searchLongestIncreasingSeq(origMatrix);
	}
}
