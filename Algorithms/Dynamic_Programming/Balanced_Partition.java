package Dynamic_Programming;

public class Balanced_Partition {
	public int balancedPartition(int array[]) {
		int minDiff = Integer.MAX_VALUE;
		int n = array.length;
		int max = 0;
		int min = Integer.MAX_VALUE;

		if (n == 0) {
			return minDiff;
		}
		if (n == 1) {
			return array[0];
		}

		for (int i = 0; i < n; i++) {
			max += array[i];
			min = (min > array[i]) ? array[i] : min;
		}

		boolean[][] auxiliary = new boolean[n][max + 1];
		for (int sum = 0; sum < max + 1; sum++) {
			if (array[0] == sum) {
				auxiliary[0][sum] = true;
			}
		}
		for (int i = 1; i < n; i++) {
			for (int sum = min; sum < max + 1; sum++) {
				if (array[i] == sum
						|| auxiliary[i - 1][sum]
						|| (sum - array[i] >= min && auxiliary[i - 1][sum
								- array[i]])) {
					auxiliary[i][sum] = true;
				}
			}
		}

		for (int division = min; division < max + 1; division++) {
			if (auxiliary[n - 1][division]) {
				if (minDiff > Math.abs(max - division - division)) {
					minDiff = Math.abs(max - division - division);
				}
			}
		}

		return minDiff;
	}

	public static void main(String args[]) {
		Balanced_Partition bp = new Balanced_Partition();
		int[] array = new int[] { 1,2,3,4,5 };
		System.out.println(bp.balancedPartition(array));
	}
}
