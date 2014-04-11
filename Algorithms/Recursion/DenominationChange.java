package Recursion;

/* 
 * 8.7 
 * Given an infinite number of quarters (25 cents), dimes (10 cents), nickels (5 cents) and pennies (1 cent), 
 * write code to calculate the number of ways of representing n cents, e.g. n = 100
 */

public class DenominationChange {
	public int makeChanges(int targetValue, int currentDenom) {
		int next_Denom = 0;
		switch (currentDenom) {
		case (25):
			next_Denom = 10;
			break;
		case (10):
			next_Denom = 5;
			break;
		case (5):
			next_Denom = 1;
			break;
		case (1):
			return 1;
		}
		int ways = 0;
		for (int i = 0; i * currentDenom <= targetValue; i++) {
			ways += makeChanges(targetValue - i * currentDenom, next_Denom);
		}
		return ways;
	}

	public static void main(String[] args) {
		DenominationChange denominationChange = new DenominationChange();
		System.out.println("# of ways = "
				+ denominationChange.makeChanges(100, 25));
	}
}
