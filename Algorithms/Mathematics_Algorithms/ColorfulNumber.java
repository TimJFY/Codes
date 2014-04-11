package Mathematics_Algorithms;

import java.util.ArrayList;
import java.util.Hashtable;

/*
 * Determine whether a number is colorful or not. 
 * 263 is a colorful number because (2,6,3,2x6,6x3,2x3x6) are all different 
 * whereas 236 is not because (2,3,6,2x3,3x6,2x3x6) have 6 twice.
 * So take all consecutive subsets of digits, 
 * take their product and ensure all the products are different.
 */

public class ColorfulNumber {

	public boolean isColorful(int num) {
		ArrayList<Integer> productSet = new ArrayList<Integer>();
		Hashtable<Integer, Integer> distinctProducts = new Hashtable<Integer, Integer>();
		int newLayerStart = 0;
		int end = 0;

		// how many products in total ? C(N+1, 2) or ∑N = (N+1)*N/2 ,
		// N is the # digital bits of num
		// compute all products from the least significant bit
		// e.g. 263: { (3), (6*3, 6), (2*6*3, 2*6, 2) }
		while (num > 0) {
			int tail = num % 10;
			num = (int) num / 10;
			end = productSet.size();
			for (int i = newLayerStart; i < end; i++) {
				int product = productSet.get(i) * tail;
				productSet.add(product);
			}
			productSet.add(tail);
			newLayerStart = end;
		}

		for (int j = 0; j < productSet.size(); j++) {
			if (distinctProducts.containsKey(productSet.get(j))) {
				return false;
			} else {
				distinctProducts.put(productSet.get(j), 1);
			}
		}
		return true;
	}

	public static void main(String args[]) {
		ColorfulNumber colorfulNumber = new ColorfulNumber();
		System.out.println(colorfulNumber.isColorful(9872));
	}
}
