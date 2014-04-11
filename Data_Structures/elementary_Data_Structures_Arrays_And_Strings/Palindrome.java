package elementary_Data_Structures_Arrays_And_Strings;

public class Palindrome {
	public int[] maxPali(char[] str) {
		int maxBegin = 0;
		int maxEnd = 0;
		int maxLen = 0;

		int center = 0;

		while (center < str.length - 1) {
			int locMaxLen = 1;
			int end = center + 1;
			int begin = center - 1;

			while (end < str.length && str[end] == str[center]) {
				end++;
				locMaxLen++;
			}
			while (begin >= 0 && end < str.length && str[begin] == str[end]) {
				begin--;
				end++;
				locMaxLen += 2;
			}
			if (locMaxLen > maxLen) {
				maxLen = locMaxLen;
				maxBegin = begin + 1;
				maxEnd = end - 1;
			}

			center++;

		}
		return new int[] { maxBegin, maxEnd, maxLen };
	}

	public static void main(String[] args) {
		String a = "sdfasdfasdfnddddn";
		char[] str = a.toCharArray();
		Palindrome palindrome = new Palindrome();
		int[] arguments = palindrome.maxPali(str);

		for (int i = arguments[0]; i <= arguments[1]; i++) {
			System.out.print(str[i] + "  ");
		}

		System.out.println();
		String b = "cbs";
		// how to determine whether a string is the substring of another one?
		int n = a.indexOf(b);
		a.contains(b);

		System.out.print("n =  " + n + "  " + a.contains(b));

	}
}
