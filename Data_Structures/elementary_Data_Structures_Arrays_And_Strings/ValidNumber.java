package elementary_Data_Structures_Arrays_And_Strings;


/*
 * http://discuss.leetcode.com/questions/241/valid-number
 */
		

public class ValidNumber {

	enum InputType {
		INVALID(0), SPACE(1), SIGN(2), DIGIT(3), DOT(4), EXPONENT(5), NUM_INPUTS(
				6);
		int value;

		private InputType(int value) {
			this.value = value;
		}
	}

	int[][] transitionMatrix = { { -1, 0, 3, 1, 2, -1 },
			{ -1, 8, -1, 1, 4, 5 }, { -1, -1, -1, 4, -1, -1 },
			{ -1, -1, -1, 1, 2, -1 }, { -1, 8, -1, 4, -1, 5 },
			{ -1, -1, 6, 7, -1, -1 }, { -1, -1, -1, 7, -1, -1 },
			{ -1, 8, -1, 7, -1, -1 }, { -1, 8, -1, -1, -1, -1 } };

	public boolean isNumber(String s) {
		int curState = 0;
		int index = 0;
		char curChar = ' ';
		while (index < s.length()) {
			int inputTypeValue = InputType.INVALID.value;
			curChar = s.charAt(index);
			if (curChar == ' ') {
				inputTypeValue = InputType.SPACE.value;
			} else if (curChar == '+' || curChar == '-') {
				inputTypeValue = InputType.SIGN.value;
			} else if (curChar >= '0' && curChar <= '9') {
				inputTypeValue = InputType.DIGIT.value;
			} else if (curChar == '.') {
				inputTypeValue = InputType.DOT.value;
			} else if (curChar == 'E' || curChar == 'e') {
				inputTypeValue = InputType.EXPONENT.value;
			}
			curState = transitionMatrix[curState][inputTypeValue];
			if (curState == -1) {
				return false;
			} else {
				index++;
			}
		}
		return curState == 1 || curState == 4 || curState == 7 || curState == 8;
	}

	public static void main(String args[]) {
		ValidNumber validNumber = new ValidNumber();
		String[] tests = { "123.05e-4", "efg123", "e33", ".98.00", ".33+44",
				"   -99e12.5", " -46456.66e1   " };
		for (String s : tests) {
			System.out.println(s + " : " + validNumber.isNumber(s));
		}

		System.out.println((char)(1+'0'));
	}
}
