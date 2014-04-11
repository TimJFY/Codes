package Recursion;

/* 
 * 8.4
 * Write a method to compute all permutations of a string(with all distinct characters)
 */
import java.util.ArrayList;

public class StringPermutation {

	public ArrayList<String> getPermutations(String tail) {
		ArrayList<String> localPermutations = new ArrayList<String>();
		if (tail == null) {
			return null;
		} else if (tail.length() == 0) {
			localPermutations.add("");
			return localPermutations;
		}

		char firstChar = tail.charAt(0);
		localPermutations = getPermutations(tail.substring(1));
		localPermutations = this.combine(localPermutations, firstChar);

		return localPermutations;
	}

	public ArrayList<String> combine(ArrayList<String> tailPermutaions,
			char singleHead) {
		ArrayList<String> constructComboPermutations = new ArrayList<String>();
		for (String permu : tailPermutaions) {
			String combo;
			// insert the head at position i
			for (int i = 0; i <= permu.length(); i++) {
				combo = permu.substring(0, i) + singleHead + permu.substring(i);
				constructComboPermutations.add(combo);
			}
		}
		return constructComboPermutations;
	}

	public void printPermutations(ArrayList<String> p) {
		int count = 1;
		for (String permu : p) {
			System.out.println(count + ": " + permu);
			count++;
		}
	}

	public static void main(String[] args) {
		String original = "abcd";

		StringPermutation stringPermutation = new StringPermutation();

		ArrayList<String> permutations = stringPermutation
				.getPermutations(original);
		stringPermutation.printPermutations(permutations);
	}
}
