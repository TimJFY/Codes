package Dynamic_Programming;

import java.util.ArrayList;

//O(n*n) version
public class Longest_Increasing_SubSequence {

	public int longset_incr_subseq(int[] sequence) {
		int n = sequence.length;
		int record[] = new int[n];
		int max = 0;
		int _max = 0;
		// extract the longest increasing subsequence
		ArrayList<Integer> index_incr_subseq = new ArrayList<Integer>();

		for (int i = 1; i < n; i++) {
			record[i] = 0;
		}
		record[0] = 1;

		for (int j = 1; j < n; j++) {
			int length_current_longest_subseq = 0;
			for (int k = 0; k < j; k++) {
				if (sequence[k] <= sequence[j]
						&& length_current_longest_subseq < record[k]) {
					length_current_longest_subseq = record[k];
				}
			}
			record[j] = length_current_longest_subseq + 1;
		}

		for (int i = 0; i < n; i++) {
			if (record[i] > max) {
				max = record[i];
			}
		}

		_max = max;
		for (int i = n - 1; i >= 0; i--) {
			if (_max == 0) {
				break;
			}
			if (record[i] == _max) {
				index_incr_subseq.add(sequence[i]);
				System.out.print(sequence[i] + " <-- ");
				_max--;
			}
		}

		return max;
	}

	public static void main(String[] args) {
		Longest_Increasing_SubSequence longest_Increasing_SubSequence = new Longest_Increasing_SubSequence();
		int[] sequence = new int[] { 2, 3, 4, 5, 6, 1, 7, 8, 4, 5, 10, 12, 24,
				78, 89, 100, 45, 32, 1, 99, 98 };
		int[] sequence2 = new int[] { 2, 3, 9, 2, 3, 5 };
		System.out.println("Length of LIS is: "
				+ longest_Increasing_SubSequence.longset_incr_subseq(sequence2));
	}
}
