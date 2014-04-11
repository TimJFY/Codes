package Mathematics_Algorithms;

import java.util.ArrayList;

public class DesirablePassword {

	ArrayList<Integer> allDesirablePasswords = new ArrayList<Integer>();

	public void findAllDesirablePasswords(int bitCount) {
		int i = 1;
		if (bitCount <= 0) {
			return;
		}
		for (int j = 0; j <= 9; j++) {
			allDesirablePasswords.add(j);
		}
		if (bitCount == 1) {
			return;
		}

		// '0' can not be the most significant bit of an integer
		allDesirablePasswords.remove(0);
		// use all N-bits passwords to create (N+1)-bits passwords, by appending
		while (i < bitCount) {
			int end = allDesirablePasswords.size();
			for (int j = 0; j < end; j++) {
				for (int tail = allDesirablePasswords.get(j) % 10 + 1; tail <= 9; tail++) {
					int newPassword = allDesirablePasswords.get(j) * 10 + tail;
					allDesirablePasswords.add(newPassword);
				}
			}
			// remove all previous passwords
			for (int k = 0; k < end; k++) {
				allDesirablePasswords.remove(0);
			}
			i++;
		}
	}

	public void showAllDesirablePasswords() {
		System.out.println("# of all desirable passwords: "
				+ allDesirablePasswords.size());
		for (Integer password : allDesirablePasswords) {
			System.out.print(password + " ");
		}
		System.out.println();
	}

	// recursive edition
	public void printAllDesirablePasswords(int number, int prev, int n) {
		if (n == 0) {
			System.out.print(number + " ");
			return;
		}
		for (int i = (prev + 1); i < (11 - n); i++) {
			printAllDesirablePasswords(number * 10 + i, i, n - 1);
		}
	}

	public static void main(String[] args) {
		DesirablePassword desirablePassword = new DesirablePassword();
		desirablePassword.findAllDesirablePasswords(3);
		desirablePassword.showAllDesirablePasswords();
		System.out.println("---------------------------------------");
		System.out.println("Recursive Method: ");
		desirablePassword.printAllDesirablePasswords(0, 0, 3);
	}

}
