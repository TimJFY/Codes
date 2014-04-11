package elementary_Data_Structures_Arrays_And_Strings;

//To determine if a string has all unique characters
public class UniqueChars {
	// for all character set whose size is no more than 256, e.g. ASCII
	public boolean isUiqueChars(String str) {
		boolean isUnique = true;
		boolean[] record = new boolean[256];

		for (int i = 0; i < str.length(); i++) {
			int loc = str.charAt(i);
			if (!record[loc]) {
				record[loc] = true;
			} else {
				isUnique = false;
				break;
			}
		}
		return isUnique;
	}

	// for all character set whose size is no more than 32, e.g. Alphabet
	public boolean isUiqueChars2(String str) {
		boolean isUnique = true;
		int currentChar = 0;
		int currentInt = 0;
		int firstChar = str.charAt(0) - 'a';
		int firstInt = 1 << firstChar;
		int accumulate = firstInt;

		for (int i = 1; i < str.length(); i++) {
			currentChar = str.charAt(i) - 'a';
			currentInt = 1 << currentChar;
			if ((accumulate & currentInt) == 0) {
				accumulate |= currentInt;
			} else {
				isUnique = false;
				break;
			}
		}
		return isUnique;
	}

	// test cases
	public static void main(String[] args) {
		String testCase1 = new String("abcdefg");
		String testCase2 = new String("abcdefga");
		String testCase3 = new String("abcdefggfsvoa");
		String testCase4 = new String("bfghjk");
		String testCase5 = new String("kljlasfao");
		String testCase6 = new String("poiouiiyuyww");
		UniqueChars uniqueChars = new UniqueChars();
		System.out
				.println("--------Method1: array(Hashing) comparison----------");
		System.out.println("testCase1 " + testCase1 + " "
				+ uniqueChars.isUiqueChars(testCase1) + "\n" + "testCase2 "
				+ testCase2 + " " + uniqueChars.isUiqueChars(testCase2) + "\n"
				+ "testCase3 " + testCase3 + " "
				+ uniqueChars.isUiqueChars(testCase3));
		System.out.println("--------Method2: bit computation----------");
		System.out.println("testCase4 " + testCase4 + " "
				+ uniqueChars.isUiqueChars2(testCase4) + "\n" + "testCase5 "
				+ testCase5 + " " + uniqueChars.isUiqueChars2(testCase5) + "\n"
				+ "testCase6 " + testCase6 + " "
				+ uniqueChars.isUiqueChars2(testCase6));
		System.out
				.println("--------Method3: sorting & (contiguous)comparison----------");
	}
}