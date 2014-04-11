package elementary_Data_Structures_Arrays_And_Strings;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

/* Implement LookAndSay function. 
 * For example, first, let user input a number, say 1. 
 * Then, the function will generate the next 10 numbers which satisfy this condition: 
 * 1, 11,21,1211,111221,312211... 
 * explanation: first number 1, second number is one 1, so 11. 
 * Third number is two 1(previous number), so 21. next number one 2 one 1, so 1211 and so on...
 */

public class NumberExtension {

	String[] allExtentions = new String[10];

	public void generateExtension(String seed) {
		String oldNumber;
		int i = 0;

		oldNumber = seed;

		while (i < 10) {
			allExtentions[i] = oldNumber;
			ArrayList<Integer> createNewNum = new ArrayList<Integer>();
			int jump = 1;
			int oneDigit = 0;
			for (int j = 0; j < oldNumber.length(); j = j + jump) {
				int ideticalIndex = 1;
				oneDigit = oldNumber.charAt(j) - '0';
				for (int k = j + 1; k < oldNumber.length(); k++) {
					if (oldNumber.charAt(k) - '0' != oneDigit) {
						break;
					}
					ideticalIndex++;
				}
				createNewNum.add(ideticalIndex);
				createNewNum.add(oneDigit);
				jump = ideticalIndex;
			}
			oldNumber = "";
			for (int l = 0; l < createNewNum.size(); l++) {
				oldNumber += createNewNum.get(l) + "";
			}
			i++;
		}
	}

	public void showExt() {
		for (int i = 0; i < 9; i++) {
			System.out.println(this.allExtentions[i]);
		}
	}

	public static void main(String[] args) {
		NumberExtension numberExtension = new NumberExtension();
		BufferedReader strin = new BufferedReader(new InputStreamReader(
				System.in));
		String seed = null;
		boolean tag = true;
		while (tag) {
			System.out.println("Input a number [0 - 9]: ");
			try {
				seed = strin.readLine();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			if ((seed.equals("")) || Integer.parseInt(seed) > 9
					|| Integer.parseInt(seed) < 0) {
				System.out.println("Wrong number: " + seed);
			} else {
				numberExtension.generateExtension(seed);
				numberExtension.showExt();
				tag = false;
			}
		}
	}
}
