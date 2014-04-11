package elementary_Data_Structures_Arrays_And_Strings;

import java.util.Hashtable;

public class RemoveRepet {
	int[] results;

	public int[] removeRepe(int[] orig) {
		Hashtable<Integer, Integer> lookupTable = new Hashtable<Integer, Integer>();
		// index for the array after removeReps
		int newArr = 0;
		for (int currentExamine = 0; currentExamine < orig.length; currentExamine++) {
			if (lookupTable.get(orig[currentExamine]) == null) {
				lookupTable.put(orig[currentExamine], 1);
				// for each word appears at first time, do in place substitution to record it 
				orig[newArr] = orig[currentExamine];
				newArr++;
			}
		}
		results = new int[newArr];
		for (int i = 0; i < newArr; i++) {
			results[i] = orig[i];
		}
		return orig;
	}

	public int[] getResult() {
		return this.results;
	}

	public static void main(String[] args) {

		RemoveRepet r = new RemoveRepet();
		int[] x = { 1, 2, 3, 4, 5, 5, 5, 5, 6, 7, 8, 8, 8, 8, 9, 9, 10, 1, 2,
				2, 3, 11, 50, 49, 6, 99 };
		int[] y;
		int[] z;
		y = r.removeRepe(x);
		z = r.getResult();
		for (int j = 0; j < y.length; j++) {
			System.out.print(y[j] + " ");
		}
		System.out.println();
		System.out
				.println("------------------------------------------------------------------------");
		for (int k = 0; k < z.length; k++) {
			System.out.print(z[k] + " ");
		}
	}
}
