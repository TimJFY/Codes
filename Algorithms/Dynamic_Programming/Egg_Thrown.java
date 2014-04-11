package Dynamic_Programming;

import java.io.BufferedReader;

import java.io.IOException;

import java.io.InputStreamReader;

public class Egg_Thrown {

	// matrix strategy[n][k] is used to store all optimal solutions to problem
	// set F(n,k)

	static int strategy[][];

	public int find_cost(int total, int egg) {

		// matrix attempts[n][k] is used to store the cost of all optimal
		// solutions to problem set F(n,k)

		int attempts[][] = new int[total + 1][egg + 1];

		strategy = new int[total + 1][egg + 1];

		for (int i = 0; i < total + 1; i++) {
			for (int j = 0; j < egg + 1; j++) {
				attempts[i][j] = Integer.MAX_VALUE;
				strategy[i][j] = 0;
			}
		}

		// F(n,1) = n
		for (int i = 1; i < total + 1; i++) {
			attempts[i][1] = i;
		}

		// F(1,n) = 1

		for (int j = 1; j < egg + 1; j++) {
			attempts[1][j] = 1;
			strategy[1][j] = 1;
		}

		for (int i = 1; i < total + 1; i++) {
			for (int j = 1; j < egg + 1; j++) {
				// r is the possible critical layer.
				// For each choice of r, there are two cases:
				// 1. if the egg is not broken, we will use k eggs to test the
				// higher (n-r) layers;
				// 2. otherwise, we should use the remaining (k-1) eggs to test
				// the lower (r-1) layers.
				// It is unknown that whether case A or case B will actually
				// happen, so we should compute the worst cost in these two
				// cases.
				// Finally, we choose the optimal strategy, which has the lowest
				// cost(times of attempts), from all choices of r.

				int criticallayer = 0;
				for (int r = 1; r < i; r++) {
					int case1 = Integer.MAX_VALUE;
					int case2 = Integer.MAX_VALUE;
					int worseOfTheTwoCases = 0;
					if (attempts[i - r][j] != Integer.MAX_VALUE)
						case1 = attempts[i - r][j] + 1;
					if (attempts[r - 1][j - 1] != Integer.MAX_VALUE)
						case2 = attempts[r - 1][j - 1] + 1;
					worseOfTheTwoCases = Math.max(case1, case2);
					if (attempts[i][j] >= worseOfTheTwoCases) {
						attempts[i][j] = worseOfTheTwoCases;
						criticallayer = r;
					}
				}

				if (criticallayer != 0)
					strategy[i][j] = criticallayer;
			}
		}

		this.show_floors(attempts, total, egg);
		return attempts[total][egg];

	}

	public void show_floors(int[][] attempts, int total, int egg) {

		int subproblem = total;
		int critical = strategy[subproblem][egg];
		int base = 0;
		int counter = 1;

		while (critical >= 1 && critical < total) {
			System.out.print("Attempt " + counter + " at layer "
					+ (critical + base) + ",");
			if (attempts[critical - 1][egg - 1] > attempts[subproblem
					- critical][egg]) {
				System.out
						.println(" if the egg is not broken, then converts into SubProblem F("
								+ (critical - 1) + "," + (egg - 1) + ") ");
				base = 0;
				subproblem = critical - 1;
				critical = strategy[critical - 1][egg - 1];
			}

			else {
				System.out
						.println(" if the egg is not broken, then converts into SubProblem F("
								+ (subproblem - critical) + "," + egg + ") ");
				base = critical + base;
				subproblem = subproblem - critical;
				critical = strategy[subproblem][egg];
			}

			System.out
					.println("-------------------------------------------------------------------------------------------");
			counter++;
		}
	}

	public static void main(String args[]) {

		int total = 0;
		int eggs = 0;
		String flag = "n";

		while (flag.equals("n")) {
			Egg_Thrown dynamic_egg = new Egg_Thrown();
			System.out.println("Please input the number of layers: ");
			
			try {
				BufferedReader strin = new BufferedReader(
						new InputStreamReader(System.in));
				total = Integer.parseInt(strin.readLine());
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			System.out.println("Please input the number of eggs: ");

			try {
				BufferedReader strin = new BufferedReader(
						new InputStreamReader(System.in));
				eggs = Integer.parseInt(strin.readLine());
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			System.out.println();

			int out = dynamic_egg.find_cost(total, eggs);

			System.out.println();
			System.out.println("Total Attempts worst case = " + out);
			System.out.println();
			System.out.println("Do you want to quit the program now? y/n");

			try {
				BufferedReader strin = new BufferedReader(
						new InputStreamReader(System.in));
				flag = strin.readLine();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			System.out.println();
		}
	}
}
