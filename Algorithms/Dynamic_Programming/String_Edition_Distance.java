package Dynamic_Programming;

public class String_Edition_Distance {

	int[][] disMatrix;
	char[][] trace;

	public int computeDistance(String orig, String target) {

		int m = orig.length();
		int n = target.length();
		disMatrix = new int[m + 1][n + 1];
		trace = new char[m + 1][n + 1];

		for (int i = 0; i <= m; i++) {
			disMatrix[i][0] = i;
			trace[i][0] = 'v';
		}
		for (int j = 0; j <= n; j++) {
			disMatrix[0][j] = j;
			trace[0][j] = 'h';
		}
		for (int i = 1; i <= m; i++) {
			for (int j = 1; j <= n; j++) {
				// no need to change, no distance
				if (orig.charAt(i - 1) == target.charAt(j - 1)) {
					disMatrix[i][j] = disMatrix[i - 1][j - 1];
					trace[i][j] = 'r';
				} else {
					// use the transformation with least cost
					int delete_append = disMatrix[i - 1][j];
					int append_delete = disMatrix[i][j - 1];
					int change = disMatrix[i - 1][j - 1];
					disMatrix[i][j] = Math.min(
							Math.min(append_delete, delete_append), change) + 1;
					if (delete_append == disMatrix[i][j] - 1) {
						trace[i][j] = 'v';
					} else if (append_delete == disMatrix[i][j] - 1) {
						trace[i][j] = 'h';
					} else {
						trace[i][j] = 'd';
					}
				}
			}
		}
		return disMatrix[m][n];
	}

	public void showSteps(String orig, String target) {
		int i = trace.length - 1;
		int j = trace[0].length - 1;
		int count = 0;
		while (i != 0 || j != 0) {
			if (trace[i][j] == 'r') {
				System.out.println("Retain '" + orig.charAt(i - 1) + "'"
						+ " At index[" + (i - 1) + "]");
				i--;
				j--;
			} else if (trace[i][j] == 'h') {
				count++;
				// append to the original string or delete from target string
				System.out.println("Append '" + target.charAt(j - 1) + "'"
						+ " to index[" + (i - 1) + "]" + "; current distance = "
						+ count);
				j--;
			} else if (trace[i][j] == 'd') {
				count++;
				// modify the original string according to the target string
				System.out
						.println("Change '" + orig.charAt(i - 1) + "'"
								+ " at index[" + (i - 1) + "]" + " into '"
								+ target.charAt(j - 1) + "'; current distance = "
								+ count);
				j--;
			} else {
				count++;
				// delete from the original string or append to target string
				System.out.println("Delete '" + orig.charAt(i - 1) + "'"
						+ " at index[" + (i - 1) + "]" + "; current distance = "
						+ count);
				i--;
			}
		}
	}

	public static void main(String[] args) {
		String_Edition_Distance string_Edition_Distance = new String_Edition_Distance();
		String a = "abcd";
		String b = "bca";
		System.out.println("Original String = '" + a + "'; Target String = '"
				+ b + "'");
		System.out.println("distance = "
				+ string_Edition_Distance.computeDistance(a, b));
		string_Edition_Distance.showSteps(a, b);
	}
}
