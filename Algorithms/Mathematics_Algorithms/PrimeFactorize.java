package Mathematics_Algorithms;

import java.util.Scanner;

public class PrimeFactorize {
	int initialFactor = 2;

	public void Factorize() {
		System.out.println("Please enter a number");
		Scanner s = new Scanner(System.in);
		int num = s.nextInt();
		System.out.println("Start Factorizing: ...");
		if (num == initialFactor)
			System.out.println(num + " = " + initialFactor);
		else {
			System.out.print(num + " = ");
			int factor = initialFactor;
			while (num > 1) {
				if (num % factor == 0) {
					num = num / factor;
					if (num > 1)
						System.out.print(factor + " * ");
					else {
						System.out.print(factor);
					}
				} else {
					factor++;
				}
			}
		}
	}

	public static void main(String[] args) {
		PrimeFactorize primeFactorize = new PrimeFactorize();
		primeFactorize.Factorize();
	}

}
