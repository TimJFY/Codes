package Recursion;

import java.util.Scanner;

/*
 * User inputs a sequence of digits. Every digit is a keystroke, that is equivalent to some character out of a sequence of characters. 
 * Digit zero and five mean NULL. The table is given below 
 0 - NULL 
 1 - v, t, f, r, q 
 2 - f, t, k 
 3 - w, z, b, g 
 4 - r, s 
 5 - NULL 
 6 - f, i, r 
 7 - p 
 8 - l, o 
 9 - p 
 * Generate all possible character sequence for a given sequence of digits. 
 * Ex - If the user input 9801, your program should generate 
 * {plv, plt, plf, plr, plq, pov, pot, pof, por, poq} (not necessarily in this order). 
 * This problem is somewhat similar to the SMS problem. I
 * t basically boils down to generating a cartesian product of the character sets corresponding to keys.
 */
public class TranslationExhaustion {

	public static String[] s = { null, "vtfrq", "ftk", "wzbg", "rs", null,
			"fir", "p", "lo", "p" };

	public void combination(String prefix, String num) {
		if (num.length() == 0)
			System.out.print(prefix + ", ");
		else {
			int x = Integer.parseInt(num.charAt(0) + "");
			String temp = new String();
			temp = s[x];
			if (temp != null) {
				for (int j = 0; j < temp.length(); j++) {
					combination(prefix + temp.charAt(j),
							num.substring(1, num.length()));
				}
			} else
				combination(prefix, num.substring(1, num.length()));
		}
	}

	public static void main(String args[]) {
		Scanner in = new Scanner(System.in);
		int num;
		System.out.print("Enter a number: ");
		num = in.nextInt();
		TranslationExhaustion translationExhaustion = new TranslationExhaustion();
		translationExhaustion.combination("", Integer.toString(num));
	}

}
