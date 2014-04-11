package Dynamic_Programming;

// !!!!!!!!!!!!!!!!!!!!NOT COMPLETED
public class Optimal_Binary_Search_Trees {

	public static double[][] expected_search_cost;
	public static double[][] search_cost_increment;
	public static int[][] roots_optimal_BST;

	public void optimal_BST(double[] probability_keys,
			double[] probability_dummies) {
		int n = probability_keys.length;
		// x-index from 1 to (n+1), y-index from 0 to n
		expected_search_cost = new double[n + 2][n + 1];
		search_cost_increment = new double[n + 2][n + 1];
		// both dimension index from 1 to n
		roots_optimal_BST = new int[n + 1][n + 1];

		for (int i = 1; i <= n + 1; i++) {
			// the subtree contains the single dummy key
			expected_search_cost[i][i - 1] = probability_dummies[i - 1];
			search_cost_increment[i][i - 1] = probability_dummies[i - 1];
		}

		for (int subtree_size = 1; subtree_size <= n; subtree_size++) {
			for (int i = 1; i <= n - subtree_size + 1; i++) {
				int j = i + subtree_size - 1;
				expected_search_cost[i][j] = Integer.MAX_VALUE;
				search_cost_increment[i][j] = search_cost_increment[i][j - 1]
						+ probability_keys[j] + probability_keys[j];
				for (int r = i; r <= j; r++) {
					double attempt_cost = expected_search_cost[i][r - 1]
							+ expected_search_cost[r + 1][j]
							+ search_cost_increment[i][j];
					if (attempt_cost > 0
							&& attempt_cost < expected_search_cost[i][j]) {
						expected_search_cost[i][j] = attempt_cost;
						roots_optimal_BST[i][j] = r;
					}
				}
			}
		}
	}

	public void construct_optimal_BST(boolean tag, int current_root, int start,
			int end) {
		if (tag) {
			System.out.println("key[" + current_root
					+ "] is the root of optimal BST");
		}
		if (current_root > start && current_root < end) {
			int left_subtree_root = roots_optimal_BST[start][current_root];
			System.out.println("key[" + left_subtree_root
					+ "] is the left child of node key[" + current_root + "]");
			construct_optimal_BST(false, left_subtree_root, start, current_root);

			int right_subtree_root = roots_optimal_BST[current_root + 1][end];
			System.out.println("key[" + right_subtree_root
					+ "] is the right child of node key[" + current_root + "]");
			construct_optimal_BST(false, right_subtree_root,
					right_subtree_root + 1, end);
		}
		if (current_root == start) {

		}
		if (current_root == end) {

		}

	}

	public static void main(String[] args) {

		double[] probability_keys = new double[] { 0.15, 0.10, 0.05, 0.10, 0.20 };
		double[] probability_dummies = new double[] { 0.05, 0.10, 0.05, 0.05,
				0.05, 0.10 };

		Optimal_Binary_Search_Trees optimal_Binary_Search_Trees = new Optimal_Binary_Search_Trees();
		optimal_Binary_Search_Trees.optimal_BST(probability_keys,
				probability_dummies);

		optimal_Binary_Search_Trees.construct_optimal_BST(true,
				roots_optimal_BST[1][roots_optimal_BST.length - 1], 1,
				roots_optimal_BST.length - 1);

	}
}
