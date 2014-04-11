package elementary_Data_Structures_Arrays_And_Strings;

// the Max sum of the subArray
public class MaximumSubarray {
	// fail to response to situation where the sequence of minsum is included in
	// the sequence of maxsum
	public int MaxSubarray_debug(int[] src) {
		int sum = src[0];
		int minsum = src[0];
		int maxsum = src[0];
		int sum1 = src[0];
		int sum2 = src[0];
		int max = src[0];

		for (int i = 1; i < src.length; i++) {
			sum += src[i];
			if (sum1 > 0) {
				sum1 += src[i];
			} else {
				sum1 = src[i];
			}

			if (sum2 < 0) {
				sum2 += src[i];
			} else {
				sum2 = src[i];
			}

			if (sum1 > maxsum) {
				maxsum = sum1;
			}

			if (sum2 < minsum) {
				minsum = sum2;
			}
		}
		max = maxsum > sum - minsum ? maxsum : sum - minsum;
		// take care if all elements are negative
		if (minsum == sum) {
			max = maxsum;
		}

		return max;
	}

	public int[] MaxSubarray(int[] src) {
		int cur_max, max;
		int i, cur_left, cur_right;
		int max_left, max_right;
		cur_max = max = src[0];
		max_left = max_right = cur_left = cur_right = 0;
		for (i = 1; i < src.length; ++i) {
			cur_max += src[i];
			if (cur_max > 0) {
				cur_right = i;
				if (max < cur_max) {
					max = cur_max;
					max_left = cur_left;
					max_right = cur_right;
				}
			} else {
				cur_max = src[i];
				// negative cases
				if (cur_max > max) {
					max = cur_max;
					max_left = cur_left = i;
					max_right = cur_right = i;
				}
				// if cur_max or src[i] is a positive value, it should be
				// recorded as the new start point
				if (cur_max < 0)
					cur_left = cur_right = i + 1;
			}
		}
		return new int[] { max, max_left, max_right };
	}

	public static void main(String[] args) {
		// int src[] = { 17, -3, 9, -2, 1, -1 };
		// int src[] = { -5, - 4, -6 };
		int src[] = { -5, -4, -3, 2, -3 };
		MaximumSubarray ms = new MaximumSubarray();
		System.out.println("The maximum sum of subarray is "
				+ ms.MaxSubarray_debug(src));
		System.out.println("------------------- ");

		System.out.println("The maximum sum of subarray is "
				+ ms.MaxSubarray(src)[0]);
		System.out.println("Index from " + ms.MaxSubarray(src)[1] + " to "
				+ ms.MaxSubarray(src)[2]);
	}

}
