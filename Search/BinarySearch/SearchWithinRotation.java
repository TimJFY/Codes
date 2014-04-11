package BinarySearch;

/*
 * 9.3
 * Given a sorted array of n integers that has been rotated an unknown number of times,give an O(logn) algorithm that finds an element in the array
 * You may assume that the array was originally sorted in increasing order
 */

public class SearchWithinRotation {

	int probe = -1;
	int targetIndex = -1;

	// to find the head of the original sorted array. O(logn)
	public void anchor(int[] rotated, int beginIndex, int endIndex) {
		int probeTemp = (beginIndex + endIndex) / 2;
		if (probeTemp > 0 && probeTemp < rotated.length) {
			if (rotated[probeTemp - 1] > rotated[probeTemp]
					&& rotated[probeTemp] < rotated[probeTemp + 1]) {
				probe = probeTemp;
			} else {
				if (rotated[beginIndex] > rotated[probeTemp]
						&& rotated[probeTemp] < rotated[endIndex]) {
					anchor(rotated, beginIndex, probeTemp);
				} else if (rotated[beginIndex] < rotated[probeTemp]
						&& rotated[probeTemp] > rotated[endIndex]) {
					anchor(rotated, probeTemp, endIndex);
				} else {
					if (rotated[beginIndex] < rotated[endIndex]) {
						probe = beginIndex;
					} else {
						probe = endIndex;
					}
				}
			}
		}
	}

	// map between the real middle and the offset middle. O(logn)
	public int binarySearchWithOffset(int[] array, int target, int offset,
			int begin, int end) {
		if (offset == -1) {
			return -1;
		}

		int realMiddle = (begin + end) / 2;
		int offsetMiddle = realMiddle + offset;

		if (offsetMiddle > array.length - 1) {
			offsetMiddle = offsetMiddle - array.length;
		}
		if (array[offsetMiddle] == target) {
			targetIndex = offsetMiddle;
		} else if (array[offsetMiddle] > target && realMiddle > 0) {
			binarySearchWithOffset(array, target, offset, begin, realMiddle);
		} else if (array[offsetMiddle] < target
				&& realMiddle < array.length - 1) {
			binarySearchWithOffset(array, target, offset, realMiddle + 1, end);
		}

		return targetIndex;
	}

	public int binarySearchWithRotation(int[] array, int target, int start,
			int end) {
		if (array.length == 0 || start > end) {
			return -1;
		}

		int mid = start + (end - start) / 2;

		if (array[mid] == target) {
			return mid;
		}
		if (array[start] <= array[mid]) {
			if (target >= array[start] && target <= array[mid]) {
				return binarySearchWithRotation(array, target, start, mid - 1);
			} else {
				return binarySearchWithRotation(array, target, mid + 1, end);
			}
		} else {
			if (target >= array[mid] && target <= array[end]) {
				return binarySearchWithRotation(array, target, mid + 1, end);
			} else {
				return binarySearchWithRotation(array, target, start, mid - 1);
			}
		}

	}

	public static void main(String[] args) {
		int[] rotated = new int[] { 15, 16, 19, 20, 25, 1, 3, 4, 5, 7, 10, 14 };
		int offset = 0;
		int target = 1;
		int target2 = 15;
		SearchWithinRotation searchWithinRotation = new SearchWithinRotation();
		searchWithinRotation.anchor(rotated, 0, rotated.length - 1);
		offset = searchWithinRotation.probe;
		System.out.println("The index of item "
				+ "'"
				+ target
				+ "'"
				+ " is "
				+ searchWithinRotation.binarySearchWithOffset(rotated, target,
						offset, 0, rotated.length - 1));
		System.out.println("----------------------------------------");
		System.out.println("The index of item "
				+ "'"
				+ target2
				+ "'"
				+ " is "
				+ searchWithinRotation.binarySearchWithRotation(rotated,
						target2, 0, rotated.length - 1));
	}
}
