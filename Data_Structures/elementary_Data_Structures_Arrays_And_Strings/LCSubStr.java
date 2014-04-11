package elementary_Data_Structures_Arrays_And_Strings;

public class LCSubStr {

	public void getLCString(char[] str1, char[] str2) {
		int i, j;
		int len1, len2;
		len1 = str1.length;
		len2 = str2.length;
		int maxLen = len1 > len2 ? len1 : len2;
		// use 3 arrays to implement the function of similarity matrix
		// e.g. a,b,c,d,e
		// a [1,0,0,0,0]
		// e [0,0,0,0,1]
		// d [0,0,0,1,0]
		// c [0,0,1,0,0]
		// d [0,0,0,2,0]
		// record the lengths of all longest sub-strings
		int[] max = new int[maxLen];
		// record the ending indexes(according to str1) of the longest sub-strings
		int[] maxIndex = new int[maxLen];
		// record the lengths of current longest sub-strings, that is, to fill
		// the similarity matrix line by line
		int[] c = new int[maxLen];

		for (i = 0; i < len2; i++) {
			for (j = len1 - 1; j >= 0; j--) {
				if (str2[i] == str1[j]) {
					// the elements in 1st row or 1st column of the matrix have
					// no precedence
					if ((i == 0) || (j == 0)) {
						c[j] = 1;
					} else
					// accumulate the length by adding 1 to its precedence(the
					// value at diagonal line), that is, c[j-1]
					{
						c[j] = c[j - 1] + 1;
					}
				} else {
					c[j] = 0;
				}
				// when find a new longest sub-string, record its length in
				// max[], its end index in maxIndex[]
				if (c[j] > max[0]) {
					max[0] = c[j];
					maxIndex[0] = j;
					// then remove all obsolete longest sub-strings
					for (int k = 1; k < maxLen; k++) {
						max[k] = 0;
						maxIndex[k] = 0;
					}
				}
				// when find another sub-string having same length with the
				// longest ones, add it
				else if (c[j] == max[0]) {
					for (int k = 1; k < maxLen; k++) {
						if (max[k] == 0) {
							max[k] = c[j];
							maxIndex[k] = j;
							break;
						}
					}

				}
			}
		}

		for (j = 0; j < maxLen; j++) {
			// output all longest sub-strings
			if (max[j] > 0) {
				System.out.print("LCSubStr " + (j + 1) + " :  ");
				for (i = maxIndex[j] - max[j] + 1; i <= maxIndex[j]; i++)
					System.out.print(str1[i]);
				System.out.println();
			}
		}

	}

	public static void main(String[] args) {
		String str1 = new String("adbba1234");
		String str2 = new String("adbbf1234sa");
		LCSubStr lcsubStr = new LCSubStr();
		lcsubStr.getLCString(str1.toCharArray(), str2.toCharArray());
	}

}