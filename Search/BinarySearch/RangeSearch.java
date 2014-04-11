package BinarySearch;

import java.util.ArrayList;

public class RangeSearch {
	static ArrayList<Integer> result = new ArrayList<Integer>();

	// O(n)
	public void rangeSearch(int[] sortedArray, int start, int end,
			int rangeStart, int rangeEnd) {
		if (start > end) {
			return;
		}
		int mid = start + (end - start) / 2;
		if (sortedArray[mid] >= rangeStart) {
			rangeSearch(sortedArray, start, mid - 1, rangeStart, rangeEnd);
		}
		if (sortedArray[mid] >= rangeStart && sortedArray[mid] <= rangeEnd) {
			result.add(sortedArray[mid]);
		}
		if (sortedArray[mid] <= rangeEnd) {
			rangeSearch(sortedArray, mid + 1, end, rangeStart, rangeEnd);
		}
	}

	// O(logn)
	public int rangeSearch_boundary(int[] sortedArray, int start, int end,
			int boundary, boolean isUpper) {

		// bound check
		if (boundary > sortedArray[sortedArray.length - 1]) {
			return sortedArray.length;
		}
		if (boundary == sortedArray[sortedArray.length - 1]) {
			if (isUpper) {
				return sortedArray.length;
			} else {
				return sortedArray.length - 2;
			}
		}
		if (boundary < sortedArray[0]) {
			return -1;
		}
		if (boundary == sortedArray[0]) {
			if (isUpper) {
				return 1;
			} else {
				return -1;
			}
		}

		// x does not exist in array, and should be placed in between
		// [Element_start, x ,Element_end] if is inserted
		if (start == end - 1) {
			if (isUpper)
				return end;
			else
				return start;
		}
		int mid = start + (end - start) / 2;
		if (sortedArray[mid] == boundary) {
			if (isUpper)
				return mid + 1;
			else
				return mid - 1;
		} else if (sortedArray[mid] > boundary) {
			return rangeSearch_boundary(sortedArray, start, mid, boundary,
					isUpper);
		} else {
			return rangeSearch_boundary(sortedArray, mid, end, boundary,
					isUpper);
		}
	}

	public static void main(String[] args) {
		RangeSearch rangeSearch = new RangeSearch();
		int[] arr = new int[] { 1, 3, 5, 7, 9, 10, 11, 18 };
		int rangeStart = 9;
		int rangeEnd = 19;
		rangeSearch.rangeSearch(arr, 0, arr.length - 1, rangeStart, rangeEnd);
		System.out.println(result);
		System.out.println("------------------------");
		result.clear();
		int lowerIndex = rangeSearch.rangeSearch_boundary(arr, 0,
				arr.length - 1, rangeStart, false);
		int upperIndex = rangeSearch.rangeSearch_boundary(arr, 0,
				arr.length - 1, rangeEnd, true);
		System.out.println("Valid indices interval: " + "[" + (lowerIndex + 1)
				+ "," + (upperIndex - 1) + "]");

		for (int i = lowerIndex + 1; i < upperIndex; i++) {
			result.add(arr[i]);

		}
		System.out.println(result);
	}
}
