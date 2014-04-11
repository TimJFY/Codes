package elementary_Data_Structures_Arrays_And_Strings;

import java.util.Hashtable;

public class Anagram {

	// all kinds of orthogonal-transformation will work, 
	// because the anagram only focuses on the quantity of each dimension, but ignore their orders
	// e.g. 1. build a 26-dimensional vector[] for each word, 
	//         and for each char(lower case), vector[char-'a']++, then compare vectors.
	// e.g  2. map each char to a distinct prime number, 
	//         char(a,b,c,d,e,...) <-> prime(2,3,5,7,11,...)
	//         compute product for each word, then compare the products.
	public boolean isAnagram(String a, String b) {
		if(a==null||b==null){
			return false;
		}
		if (a.length() != b.length()) {
			return false;
		}
		// arrange the information of string a
		char[] CharArr_a = a.toCharArray();
		int distinct_elements_a = 0;
		Hashtable<Character, Integer> lookupTable = new Hashtable<Character, Integer>();
		for (int i = 0; i < a.length(); i++) {
			if (!lookupTable.containsKey(CharArr_a[i])) {
				lookupTable.put(CharArr_a[i], 1);
				distinct_elements_a++;
			} else {
				lookupTable
						.put(CharArr_a[i], lookupTable.get(CharArr_a[i]) + 1);
			}
		}
		// match with b
		char[] CharArr_b = b.toCharArray();
		int matched_elements = 0;
		for (int j = 0; j < b.length(); j++) {
			// new element in b or already matched element in a
			if (!lookupTable.containsKey(CharArr_b[j])
					|| lookupTable.get(CharArr_b[j]) == 0) {
				return false;
			} else {
				lookupTable
						.put(CharArr_b[j], lookupTable.get(CharArr_b[j]) - 1);
				if (lookupTable.get(CharArr_b[j]) == 0) {
					matched_elements++;
					if (matched_elements == distinct_elements_a) {
						return j == b.length() - 1;
					}
				}
			}
		}
		return false;
	}
	
	public static void main(String[] args){
		Anagram anagram=new Anagram();
		String a="abccccde";
		String b="cccedbac";
		System.out.println(anagram.isAnagram(a, b));
	}
}
