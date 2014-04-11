package Dynamic_Programming;

public class MatrixChain_Mutiplication {

	public static int cost[][];
	public static int index[][];

	public void matrix_chain_order(int[] demsion_seq) {
		int n = demsion_seq.length;
		cost = new int[n][n];
		index = new int[n][n];

		for (int diagnal = 1; diagnal < n; diagnal++) {
			cost[diagnal][diagnal] = 0;
		}
		// l is the chain length
		for (int l = 2; l < n; l++) {
			for (int i = 1; i < n - l + 1; i++) {
				int j = i + l - 1;
				cost[i][j] = Integer.MAX_VALUE;
				// test all (j-i) possible combinations
				for (int k = i; k < j; k++) {
					int current_cost = cost[i][k] + cost[k + 1][j]
							+ demsion_seq[i - 1] * demsion_seq[k]
							* demsion_seq[j];
					if (current_cost >= 0 && current_cost < cost[i][j]) {
						cost[i][j] = current_cost;
						index[i][j] = k;
					}
				}
			}
		}
	}

	public void print_optiamal_parens(int[][] index, int i, int j) {
		if (i == j) {
			System.out.print("A" + i);
		} else {
			System.out.print("(");
			print_optiamal_parens(index, i, index[i][j]);
			print_optiamal_parens(index, index[i][j] + 1, j);
			System.out.print(")");
		}
	}

	public static void main(String[] args) {
		int[] demsion_seq = new int[] { 30, 35, 15, 5, 10, 20, 25 };
		MatrixChain_Mutiplication matrixChain_Mutiplication = new MatrixChain_Mutiplication();
		matrixChain_Mutiplication.matrix_chain_order(demsion_seq);
		matrixChain_Mutiplication.print_optiamal_parens(index, 1,
				demsion_seq.length - 1);
		//Optimal subStructure
		System.out.println();
		System.out.println("--------------------");
		matrixChain_Mutiplication.print_optiamal_parens(index, 2 ,demsion_seq.length - 2);
	}
}
