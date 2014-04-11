package Dynamic_Programming;

import java.util.ArrayList;
import java.util.StringTokenizer;

//print neatly problem
public class Print_Neatly {

	public static int total_Space[];
	public static int start_Index[];
	public static int estimated_Space[][];
	public static int line_Space[][];
	// if the length of a word is larger than the line size, printing will fail.
	static int wordTooLong = -1;

	public int[] print_satisitics(ArrayList<String> text, int lineSize) {
		int num_Words = text.size();
		total_Space = new int[num_Words];
		estimated_Space = new int[num_Words][num_Words];
		start_Index = new int[num_Words];
		line_Space = new int[num_Words][num_Words];

		for (int i = 1; i < num_Words; i++) {
			estimated_Space[i][i] = lineSize - text.get(i).length();
			if (estimated_Space[i][i] < 0) {
				wordTooLong = i;
			}
			for (int j = i + 1; j < num_Words; j++) {
				estimated_Space[i][j] = estimated_Space[i][j - 1]
						- text.get(j).length() - 1;
			}
		}

		for (int i = 1; i < num_Words; i++) {
			for (int j = i; j < num_Words; j++) {
				if (estimated_Space[i][j] < 0) {
					line_Space[i][j] = Integer.MAX_VALUE;
				} else if (j == num_Words) {
					line_Space[i][j] = 0;
				} else {
					line_Space[i][j] = (int) Math.pow(estimated_Space[i][j], 3);
				}
			}
		}

		total_Space[0] = 0;
		for (int j = 1; j < num_Words; j++) {
			total_Space[j] = Integer.MAX_VALUE;
			for (int i = 1; i <= j; i++) {
				if (line_Space[i][j] != Integer.MAX_VALUE
						&& total_Space[j] > total_Space[i - 1]
								+ line_Space[i][j]) {
					total_Space[j] = total_Space[i - 1] + line_Space[i][j];
					start_Index[j] = i;
				}
			}
		}

		return start_Index;
	}

	public void print_out(ArrayList<String> text, int lineSize,
			int[] start_Index, int end_Index_of_One_Line) {
		int start_Index_of_One_Line = start_Index[end_Index_of_One_Line];
		if (start_Index_of_One_Line > 1) {
			print_out(text, lineSize, start_Index, start_Index_of_One_Line - 1);
		}
		print_With_Space(text, lineSize, start_Index_of_One_Line,
				end_Index_of_One_Line);
	}

	public void print_With_Space(ArrayList<String> text, int lineSize,
			int start, int end) {
		int used_Space = 0;
		for (int i = start; i <= end; i++) {
			if (text.get(i).length() + used_Space < lineSize) {
				System.out.print(text.get(i) + "*");
			} else {
				System.out.print(text.get(i));
			}
			used_Space = used_Space + text.get(i).length() + 1;
		}
		for (int j = 1; j <= lineSize - used_Space; j++) {
			System.out.print("*");
		}
		System.out.println();
	}

	public static void main(String[] args) {
		String text = "This is a test for print neatly problem using dynamic programming Its original description is the problem 4 in chapter 15 of INTRODUCTION TO ALGORITHMS 3rd Edition";
		int lineSize = 15;
		// int lineSize = 10;

		ArrayList<String> test = new ArrayList<String>();
		StringTokenizer tokens = new StringTokenizer(text, " ");
		test.add("BASE");
		while (tokens.hasMoreTokens()) {
			test.add(tokens.nextToken());
		}

		Print_Neatly print_Neatly = new Print_Neatly();
		print_Neatly.print_satisitics(test, lineSize);
		if (wordTooLong != -1) {
			System.out
					.println("The "
							+ wordTooLong
							+ "th word is too long to fit in a single line, Printing can not be carried out.");
		} else {
			System.out.println("Now Printing: (each * represents a space)");
			System.out.println("------------------------------------- ");
			print_Neatly
					.print_out(test, lineSize, start_Index, test.size() - 1);
		}
	}

}
