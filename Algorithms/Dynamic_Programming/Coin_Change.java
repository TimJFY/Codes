package Dynamic_Programming;

import java.util.ArrayList;

/*
 * Coin change is the problem of finding the number of ways to make change for a target amount given a set of denominations.  
 * It is assumed that there is an unlimited supply of coins for each denomination.  
 * An example will be finding change for target amount 4 using change of 1,2,3 for which the solutions are (1,1,1,1), (2,2), (1,1,2), (1,3).  
 * As you can see, the optimal solution can be (2,2) or (1,3).  
 * The attached Java program solves both the problems of "find all combinations" and "find the optimal solution (which takes the least number of coins)".
 */

public class Coin_Change {
	// contains the optimal solution during every recurrence step.
	public int OPT[][];
	// string representation of optimal solutions.
	public String optimalChange[][];
	// List of all possible solutions for the target
	public ArrayList<String> allPossibleChanges = new ArrayList<String>();
	// The target amount.
	private int target;
	// Copy of the denominations that was used to generate this solution
	public int[] denoms;

	public Coin_Change(int target, int[] denoms) {
		this.target = target;
		this.denoms = denoms;
		OPT = new int[denoms.length][target + 1];
		optimalChange = new String[denoms.length][target + 1];

	}

	/**
	 * Find all possible solutions recursively
	 * 
	 * @param tsoln
	 *            - The current solution string
	 * @param startIx
	 *            - The start index in the denomination array.
	 * @param remainingTarget
	 *            - The remaining target amount to be satisfied.
	 */
	private void findAllCombinationsRecursive(String tsoln, int startIx,
			int remainingTarget) {
		for (int i = startIx; i < this.denoms.length; i++) {
			int temp = remainingTarget - this.denoms[i];
			String tempSoln = tsoln + "" + this.denoms[i] + ",";
			if (temp < 0) {
				break;
			}
			if (temp == 0) {
				// reached the answer hence quit from the loop
				this.allPossibleChanges.add(tempSoln);
				break;
			} else {
				// target not reached, try the solution recursively with the
				// current denomination as the start point.
				findAllCombinationsRecursive(tempSoln, i, temp);
			}
		}
	}

	/**
	 * Find the optimal solution for a given target value and the set of
	 * denominations
	 */
	public void findOptimalChange() {

		StringBuilder sb = new StringBuilder();

		// initialize the solution structure
		for (int i = 0; i < this.OPT[0].length; i++) {
			this.OPT[0][i] = i;
			this.optimalChange[0][i] = sb.toString();
			sb.append(denoms[0] + " ");
		}

		// Read through the following for more details on the explanation
		// of the algorithm.
		// http://condor.depaul.edu/~rjohnson/algorithm/coins.pdf
		for (int i = 1; i < denoms.length; i++) {
			for (int j = 0; j < target + 1; j++) {
				int value = j;
				int targetWithPrevDenomiation = this.OPT[i - 1][j];
				int ix = (value) - denoms[i];
				if (ix >= 0 && (denoms[i] <= value)) {
					int x2 = denoms[i] + this.OPT[i][ix];
					if (x2 <= target
							&& (1 + this.OPT[i][ix] < targetWithPrevDenomiation)) {
						String temp = this.optimalChange[i][ix] + denoms[i]
								+ " ";
						this.optimalChange[i][j] = temp;
						this.OPT[i][j] = 1 + this.OPT[i][ix];
					} else {
						this.optimalChange[i][j] = this.optimalChange[i - 1][j]
								+ " ";
						this.OPT[i][j] = targetWithPrevDenomiation;
					}
				} else {
					this.optimalChange[i][j] = this.optimalChange[i - 1][j];
					this.OPT[i][j] = targetWithPrevDenomiation;
				}
			}
		}
	}

	public void showResults() {
		System.out.println("The demominations are: ");
		for (int i = 0; i < denoms.length; i++) {
			System.out.print(denoms[i] + " ");
		}
		System.out.println();
		System.out.println("The target amount is: ");
		System.out.println(target);
		System.out.println("-------------------------------------");
		System.out.println("All change strategies are as follows: ");
		for (int i = 0; i < allPossibleChanges.size(); i++) {
			System.out.println((i + 1) + " : " + allPossibleChanges.get(i));
		}
		System.out.println("The optimal change strategy is: ");
		System.out.println(optimalChange[denoms.length - 1][target]);
	}

	public static void main(String[] args) {
		int target = 12;
		int[] denoms = new int[] { 1, 6, 10 };

		Coin_Change soln = new Coin_Change(target, denoms);
		String tempSoln = new String();
		soln.findAllCombinationsRecursive(tempSoln, 0, target);
		soln.findOptimalChange();
		soln.showResults();
	}
}
