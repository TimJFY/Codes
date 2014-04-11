package Recursion;

/* 
 * 8.5
 * Implement an algorithm to print all valid (e g , properly opened and closed) combinations of n-pairs of parentheses
 */

import java.util.ArrayList;

public class ParenthesesPairs {
	public boolean printParenthesesPairs(ArrayList<Character> seq, int index,
			int leftRemain, int rightRemain) {
		boolean tag = false;
		if (leftRemain < 0 || rightRemain < leftRemain) {
			return false;
		}
		if (leftRemain == 0 && rightRemain == 0) {
			System.out.println(seq);
			return true;
		}
		if (leftRemain > 0) {
			seq.set(index, '{');
			tag = printParenthesesPairs(seq, index + 1, leftRemain - 1,
					rightRemain);
		}
		if (rightRemain > leftRemain) {
			seq.set(index, '}');
			tag = printParenthesesPairs(seq, index + 1, leftRemain,
					rightRemain - 1);
		}

		return tag;
	}

	public static void main(String[] args) {
		int pairs = 3;
		ArrayList<Character> seq = new ArrayList<Character>();
		for (int i = 0; i < pairs * 2; i++) {
			seq.add('|');
		}
		ParenthesesPairs parenthesesPairs = new ParenthesesPairs();
		parenthesesPairs.printParenthesesPairs(seq, 0, pairs, pairs);
	}
}
