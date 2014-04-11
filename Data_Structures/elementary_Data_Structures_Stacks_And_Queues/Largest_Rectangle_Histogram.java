package elementary_Data_Structures_Stacks_And_Queues;

import java.util.Stack;

public class Largest_Rectangle_Histogram {
	public int largestRectangleArea(int[] height) {

		Stack<Integer> heightIndices = new Stack<Integer>();
		int i = 0;
		int result = 0;
		while (i < height.length) {
			while (heightIndices.size() > 0
					&& height[heightIndices.peek()] >= height[i]) {
				int j = heightIndices.pop();
				int curMax = 0;
				if (heightIndices.size() > 0) {
					curMax = (i - heightIndices.peek() - 1) * height[j];
				} else {
					curMax = i * height[j];
				}
				result = result > curMax ? result : curMax;
			}
			heightIndices.push(i);
			i++;
		}

		while (heightIndices.size() > 0) {
			int j = heightIndices.pop();
			int curMax = 0;
			if (heightIndices.size() > 0) {
				curMax = (i - heightIndices.peek() - 1) * height[j];
			} else {
				curMax = i * height[j];
			}
			result = result > curMax ? result : curMax;
		}

		return result;
	}

	public static void main(String[] args) {
		Largest_Rectangle_Histogram lrh = new Largest_Rectangle_Histogram();
		int[] height = new int[] { 0, 0 };
		System.out.println(lrh.largestRectangleArea(height));
	}
}
