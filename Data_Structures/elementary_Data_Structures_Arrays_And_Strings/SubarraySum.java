package elementary_Data_Structures_Arrays_And_Strings;

/*
 * Substring Addition 
 * Write a program to add the substring 
 * eg :say you have a list {1 7 6 3 5 8 9 } and user is entering a sum 16.
 * Output should display (2-4),(3-5) that is {7, 6, 3},{3, 5, 8} cause 7+6+3=3+5+8=16.
 */
public class SubarraySum {

	// O(n2)): use the  
	public void SubArraySum(int sum, int[] input) {
		int[] suffixes = new int[input.length];
		int curSum = 0;
		for (int i = input.length - 1; i >= 0; i--) {
			curSum += input[i];
			suffixes[i] = curSum;
		}
		for (int i = 0; i < suffixes.length; i++) {
			for (int j = i; j < suffixes.length; j++) {
				if (suffixes[i] - suffixes[j] == sum) {
					print(i, j, input);
				}
			}
		}
	}

	public void print(int start, int end, int[] input) {
		for (int k = start; k < end; k++) {
			System.out.print(input[k] + ",");
		}
		System.out.println();
	}

	public static void main(String[] args) {
		SubarraySum subarraySum = new SubarraySum();
		int[] input = new int[] { 1, 7, 6, 3, 5, 8, 9 };
		int sum = 16;
		subarraySum.SubArraySum(sum, input);
	}
}
