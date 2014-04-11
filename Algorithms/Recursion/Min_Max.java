package Recursion;

public class Min_Max {

	public int[] getMinMax(int[] array, int start, int end) {

		int result[] = new int[] { Integer.MAX_VALUE, Integer.MIN_VALUE };
		int resultFir[] = new int[2];
		int resultSec[] = new int[2];

		resultFir[0] = array[start];
		resultFir[1] = array[start];
		resultSec[0] = array[start];
		resultSec[1] = array[start];

		if (end - start > 0) {
			resultFir = getMinMax(array, start, start + (end - start) / 2);
			resultSec = getMinMax(array, start + (end - start) / 2 + 1, end);
		}
		if (result[0] > resultFir[0])
			result[0] = resultFir[0];
		if (result[0] > resultSec[0])
			result[0] = resultSec[0];
		if (result[1] < resultFir[1])
			result[1] = resultFir[1];
		if (result[1] < resultSec[1])
			result[1] = resultSec[1];

		return result;
	}

	public static void main(String[] args) {
		int[] test1 = new int[] { -1, -9, 88, 3, 4, 7, 2, 3, 33, 2, 0, 55, 91,
				27, 10000, -90 };
		//int[] test = new int[] { 1, 9, 10, 11 };
		Min_Max min_Max = new Min_Max();
		int[] s = min_Max.getMinMax(test1, 0, test1.length - 1);
		System.out.println("Min = "+s[0] + " , Max= " + s[1]);
	}
}
