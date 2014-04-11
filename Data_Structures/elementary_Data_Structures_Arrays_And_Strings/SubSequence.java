package elementary_Data_Structures_Arrays_And_Strings;

public class SubSequence {

	public boolean isSubSeq(String a, String b) {

		if (a == null) {
			return true;
		} else if (b == null) {
			return false;
		}

		if (a.length() > b.length()) {
			return false;
		}

		int jump = 0;
		for (int i = 0; i < a.length(); i++) {
			for (int j = jump;; j++) {
				// when find a matching in string b, the next matching search should
				// start from the currently matched char
				if (a.charAt(i) == b.charAt(j)) {
					jump = j + 1;
					if (j == b.length() - 1) {
						// reaching the ends of both a and b
						if (i == a.length() - 1) {
							return true;
						}
						// b is used up, but a is not completely matched
						return false;
					}
					break;
				}
				// can not find any matching char in b
				if (j == b.length() - 1) {
					return false;
				}
			}
		}
		return true;
	}

	public static void main(String[] args) {
		SubSequence subSequence = new SubSequence();
		System.out.println(subSequence.isSubSeq("xaxbxb", "xxxxxbaxxxbxa"));
	}
}
