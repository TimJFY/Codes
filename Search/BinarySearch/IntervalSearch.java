package BinarySearch;

import java.util.ArrayList;

public class IntervalSearch {
	public static ArrayList<Integer> results = new ArrayList<Integer>();
	static int minIndex = -1;
	static int maxIndex = -1;

	public void intervalSerach(int[] array, int target, int start, int end) {
		if (start <= end) {
			int mid = start + (end - start) / 2;
			if (array[mid] == target) {
				intervalSerach(array, target, start, mid - 1);
				intervalSerach(array, target, mid + 1, end);
				minIndex = (minIndex > mid || minIndex == -1) ? mid : minIndex;
				maxIndex = (maxIndex < mid || maxIndex == -1) ? mid : maxIndex;
				results.add(mid);
			} else if (array[mid] > target) {
				intervalSerach(array, target, start, mid - 1);
			} else {
				intervalSerach(array, target, mid + 1, end);
			}
		}
	}

	public static void main(String args[]) {
		IntervalSearch intervalSearch = new IntervalSearch();
		int[] arr = new int[] { 1,2,2,2,2,2,2,2,6 };
		intervalSearch.intervalSerach(arr, 2, 0, arr.length - 1);
		System.out.println(results);
		System.out.println("Valid Interval: " + "[" + minIndex + "," + maxIndex
				+ "]");
	}
}
