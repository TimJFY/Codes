package Mathematics_Algorithms;

import java.util.Set;

public class NextPermutation {

	long result;
	boolean isIncrease;
	long pivotDigitOffest;

	public NextPermutation() {
		result = Long.MAX_VALUE;
		isIncrease = true;
		pivotDigitOffest = 0;
	}

	// Step 1: find the minimal number that is larger than the original one
	// only by trying to switch 2 digits
	public void nextPermutation(long num) {
		result = Long.MAX_VALUE;
		isIncrease = true;
		pivotDigitOffest = 0;
		long digitProbe = 1;
		long medianResult = Long.MAX_VALUE;
		while (digitProbe < num) {
			long compareDigProbe = digitProbe * 10;
			// value of current digit
			long curDigit = (num % compareDigProbe) / digitProbe;
			while (compareDigProbe < num) {
				long preDigit = (num % (compareDigProbe * 10))
						/ compareDigProbe;
				if (curDigit <= preDigit) {
					compareDigProbe *= 10;
				} else {
					if (swapDigits(curDigit, digitProbe, preDigit,
							compareDigProbe, num) < medianResult) {
						medianResult = swapDigits(curDigit, digitProbe,
								preDigit, compareDigProbe, num);
						isIncrease = false;
						pivotDigitOffest = compareDigProbe;
					}
					// for current digit, its nearest switch has been completed,
					// go next(higher) digit, The whole loop CANNOT be
					// terminated here.
					break;
				}
			}
			digitProbe *= 10;
		}

		if (!isIncrease) {
			result = adjustOrder(medianResult, pivotDigitOffest);
		} else {
			// the original number is already the largest in the permutation, so
			// return the smallest one. you can also implement it by reversing
			// the original number
			result = adjustOrder(num, digitProbe);
		}
	}

	public long swapDigits(long digit1, long offset1, long digit2,
			long offset2, long num) {
		return num + (digit1 - digit2) * (offset2 - offset1);
	}

	// Step 2: adjust all digits after the pivot(the remainder) into increasing
	// order
	public long adjustOrder(long num, long pivotDigitOffest) {
		long remainder = num % pivotDigitOffest;
		num -= remainder;

		while (true) {
			long newRemainder = remainder;
			long digitProbe = 1;
			long compareDigProbe = digitProbe * 10;
			long curDigit = (newRemainder % compareDigProbe) / digitProbe;
			// for each digit, keep on switching if capable
			while (compareDigProbe < newRemainder) {
				long preDigit = (newRemainder % (compareDigProbe * 10))
						/ compareDigProbe;
				if (curDigit < preDigit) {
					newRemainder = swapDigits(curDigit, digitProbe, preDigit,
							compareDigProbe, newRemainder);
					// since the current digit has been switched to a new
					// position, update its offset
					digitProbe = compareDigProbe;
				}
				compareDigProbe *= 10;
			}
			// can not generate smaller remainder <--> increasing order is built
			if (newRemainder == remainder) {
				break;
			} else {
				remainder = newRemainder;
			}
		}

		return num + remainder;
	}

	// public void nextPermutation(int[] num) {
	//
	// if (num.length != 0) {
	//
	// int[] medianNum = new int[num.length];
	// int[] temp = num.clone();
	// for (int i = 0; i < medianNum.length; i++) {
	// medianNum[i] = Integer.MAX_VALUE;
	// }
	//
	// boolean isIncrease = true;
	// int pivotDigit = 0;
	// int i = num.length - 1;
	// while (i >= 0) {
	// int j = i - 1;
	// while (j >= 0) {
	// temp = num.clone();
	// if (num[j] < num[i]) {
	// int switchTemp = temp[i];
	// temp[i] = temp[j];
	// temp[j] = switchTemp;
	// isIncrease = false;
	// pivotDigit = j + 1;
	// if (arrayCompare(medianNum, temp)) {
	// medianNum = temp.clone();
	// }
	// break;
	// }
	// j--;
	// }
	// i--;
	// }
	// if (isIncrease) {
	// medianNum = num.clone();
	// }
	// num = adjustOrder(medianNum, pivotDigit);
	// System.out.println(num[0]);
	// }
	// }
	//
	// public boolean arrayCompare(int[] a, int[] b) {
	// boolean result = false;
	// for (int i = 0; i < a.length; i++) {
	// if (a[i] > b[i]) {
	// return true;
	// } else if (a[i] < b[i]) {
	// return false;
	// }
	// }
	// return result;
	// }
	//
	// public int[] adjustOrder(int[] num, int startDigit) {
	// int i = startDigit;
	// while (i < num.length) {
	// int curMin = i;
	// int j = i + 1;
	// while (j < num.length) {
	// if (num[j] < num[curMin]) {
	// curMin = j;
	// }
	// j++;
	// }
	// int switchTemp = num[curMin];
	// num[curMin] = num[i];
	// num[i] = switchTemp;
	// i++;
	// }
	// return num;
	// }

	public void nextPermutation(int[] num) {
		if (num.length != 0) {
			int targetDigitPosition = -1;
			int originDigit = 0;
			int originDigitPostion = 0;

			boolean isIncrease = true;
			int i = num.length - 1;
			while (i >= 0) {
				int j = i - 1;
				while (j >= 0) {
					if (num[j] < num[i]) {
						isIncrease = false;
						if (j > targetDigitPosition) {
							targetDigitPosition = j;
							originDigit = num[i];
							originDigitPostion = i;
						} else if (j == targetDigitPosition
								&& originDigit > num[i]) {
							originDigit = num[i];
							originDigitPostion = i;
						}
						break;
					}
					j--;
				}
				i--;
			}

			if (!isIncrease) {
				int temp = num[originDigitPostion];
				num[originDigitPostion] = num[targetDigitPosition];
				num[targetDigitPosition] = temp;
			}
			num = adjustOrder(num, targetDigitPosition + 1);
			// System.out.println(num);
		}
	}

	public int[] adjustOrder(int[] num, int startDigit) {
		int i = startDigit;
		while (i < num.length) {
			int curMin = i;
			int j = i + 1;
			while (j < num.length) {
				if (num[j] < num[curMin]) {
					curMin = j;
				}
				j++;
			}
			int switchTemp = num[curMin];
			num[curMin] = num[i];
			num[i] = switchTemp;
			i++;
		}
		return num;
	}

	
	
	public static void main(String args[]) {
		NextPermutation nextPermutation = new NextPermutation();

		nextPermutation.nextPermutation(1111);
		System.out.println(1111 + " ---> " + nextPermutation.result);
		nextPermutation.nextPermutation(431);
		System.out.println(431 + " ---> " + nextPermutation.result);
		nextPermutation.nextPermutation(4679);
		System.out.println(4679 + " ---> " + nextPermutation.result);
		nextPermutation.nextPermutation(1895);
		System.out.println(34326778 + " ---> " + nextPermutation.result);
		nextPermutation.nextPermutation(944984326);
		System.out.println(944984326 + " ---> " + nextPermutation.result);
		nextPermutation.nextPermutation(17822943498765L);
		System.out.println(17822943498765L + " ---> " + nextPermutation.result);

		nextPermutation.nextPermutation(new int[] { 1, 3, 2 });
		
		

	}
	
}
