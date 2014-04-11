package Mathematics_Algorithms;

public class CountArithSeq {
	public int solution(int[] A) {
		// write your code in Java SE 7
		// length of the array
		int n = A.length;
		int result = 0;
		// number of maximum results
		int boundary = 1000000000;
		// impossible to form any arithmetic slice
		if (n <= 2) {
			return result;
		}
		// In this array, element at index of i stores the number of arithmetic
		// slices that are ended at ith element in the original array
		int[] seqCount = new int[n];
		// each time, we examine 3 consecutive elements
		int start = 0;
		int end = start + 2;
		// last but 2
		int pivot = end - 1;
		while (end < n) {
			long first = A[start];
			long last = A[end];
			long middle = (long) A[pivot] * 2;
			// can form an arithmetic slices
			if (middle == first + last) {
				seqCount[end] = seqCount[end - 1] + 1;
				end++;
				start++;
				pivot++;
			}
			// if not, update start
			else {
				start = end - 1;
				end = end + 1;
				pivot = end - 1;
			}
		}

		// sum up all possibilities
		for (int i = 0; i < n; i++) {
			result += seqCount[i];
			if (result > boundary) {
				return -1;
			}
		}
		return result;
	}

	public static void main(String args[]) {
		CountArithSeq cas = new CountArithSeq();
		System.out.println(cas.solution(new int[] { 1, 1, 1, 2, 3, 4,
				-2147483648, -2147483648, -2147483648 }));
	}
}
