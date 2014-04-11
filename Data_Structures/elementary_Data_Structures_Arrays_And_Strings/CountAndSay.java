package elementary_Data_Structures_Arrays_And_Strings;

public class CountAndSay {
	public String countAndSay(int n) {
		String curSeq = "1";
		String nexSeq = "";
		if (n == 1) {
			return curSeq;
		}
		for (int i = 1; i < n; i++) {
			char curDigit = curSeq.charAt(0);
			int index = 0;
			while (index < curSeq.length()) {
				int count = 0;
				while (index < curSeq.length()
						&& curSeq.charAt(index) == curDigit) {
					count++;
					index++;
				}
				nexSeq += (char) (count + '0');
				nexSeq += curDigit;
				if (index < curSeq.length()) {
					curDigit = curSeq.charAt(index);
				}
			}
			curSeq = nexSeq;
			nexSeq = "";
		}
		return curSeq;
	}
}
