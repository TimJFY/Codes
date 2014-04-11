package Dynamic_Programming;

//longest common sequence - memorization version
public class Longest_Common_Sequence_Memo {

	public void lcs_length_memo(char[] X, char[] Y) {
		int m = X.length;
		int n = Y.length;
		int look_up[][] = new int[m][n];
		for (int i = 1; i < m; i++) {
			look_up[i][0] = 0;
		}
		for (int j = 0; j < n; j++) {
			look_up[0][j] = 0;
		}
		for (int i = 1; i < m; i++) {
			for (int j = 1; j < n; j++) {
				look_up[i][j] = -1;
			}
		}

		// the first element is a padding, ignore it
		look_up_lcs(X, Y, look_up, m - 1, n - 1);
		
		print_lcs(X, Y, look_up, m - 1, n - 1);
	}

	public int look_up_lcs(char[] X, char[] Y, int[][] look_up, int i, int j) {
		int lcs_len;
		if (look_up[i][j] != -1) {
			return look_up[i][j];
		}
		if (i == 0 || j == 0) {
			lcs_len = 0;
		} else if (X[i] == Y[j]) {
			lcs_len = look_up_lcs(X, Y, look_up, i - 1, j - 1) + 1;
		} else {
			lcs_len = Math.max(look_up_lcs(X, Y, look_up, i, j - 1),
					look_up_lcs(X, Y, look_up, i - 1, j));
		}

		look_up[i][j] = lcs_len;

		return look_up[i][j];
	}

	public void print_lcs(char[] X, char[] Y, int[][] look_up, int i, int j) {
		if (i == 0 || j == 0) {
			System.out.println("Reach the beginning.");
		} else {
			if (X[i] == Y[j]) {
				print_lcs(X, Y, look_up, i - 1, j - 1);
				System.out.print(X[i]+"---->");
			} else if (look_up[i - 1][j] >= look_up[i][j - 1]) {
				print_lcs(X, Y, look_up, i - 1, j);
			} else {
				print_lcs(X, Y, look_up, i, j - 1);
			}
		}
	}

	public static void main(String[] args) {
		Longest_Common_Sequence_Memo look_Up_LCS_Memo = new Longest_Common_Sequence_Memo();
		// * - padding
		char[] X = new char[] { '*', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h',
				'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q' };
		char[] Y = new char[] { '*', 'a', 'c', 'b', 'd', 'e', 'g', 'i', 'f',
				'i', 'j', 'k', 'm', 'n', 'l', 'q', 'o', 'p' };
		look_Up_LCS_Memo.lcs_length_memo(X, Y);
	}
}
