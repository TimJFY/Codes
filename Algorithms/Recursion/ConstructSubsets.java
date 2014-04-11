package Recursion;

import java.util.ArrayList;

/* 
 * 8.3 
 * Write a method that returns all subsets of a set
 */
public class ConstructSubsets {

	public ArrayList<ArrayList<Integer>> findAllSubsets(ArrayList<Integer> set,
			int index) {
		ArrayList<ArrayList<Integer>> allSubSets = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Integer>> allSubSetsExtented = new ArrayList<ArrayList<Integer>>();
		if (index == set.size()) {
			ArrayList<Integer> emptySet = new ArrayList<Integer>();
			allSubSets.add(emptySet);
			return allSubSets;
		}
		Integer firstItem = set.get(index);
		allSubSets = findAllSubsets(set, index + 1);
		for (ArrayList<Integer> subSet : allSubSets) {
			// Alternative
			// ArrayList<Integer> copy = (ArrayList<Integer>) subSet.clone();
			ArrayList<Integer> copy = new ArrayList<Integer>();
			copy.addAll(subSet);
			copy.add(firstItem);
			allSubSetsExtented.add(copy);
		}
		allSubSets.addAll(allSubSetsExtented);
		return allSubSets;
	}

	// most significant bit corresponds to the last element e.g. 001 ---> [, ,3]
	public ArrayList<ArrayList<Integer>> findAllSubsetsNoRecursion(
			ArrayList<Integer> set) {
		ArrayList<ArrayList<Integer>> allSubSets = new ArrayList<ArrayList<Integer>>();
		int bitDiagram = (1 << set.size()) - 1;
		for (int bitValue = bitDiagram; bitValue >= 0; bitValue--) {
			ArrayList<Integer> oneSubSet = new ArrayList<Integer>();
			int bitPointer = 0;
			for (int j = bitValue; j != 0; j = j >> 1) {
				if ((j & 1) == 1) {
					oneSubSet.add(set.get(bitPointer));
				}
				bitPointer++;
			}
			allSubSets.add(oneSubSet);
		}
		return allSubSets;
	}

	public static void main(String[] args) {
		ArrayList<Integer> set = new ArrayList<Integer>();
		set.add(1);
		set.add(2);
		set.add(3);
		set.add(4);
		//set.add(5);
		int count = 1;
		ConstructSubsets constructSubsets = new ConstructSubsets();
		 ArrayList<ArrayList<Integer>> allSubSets = constructSubsets.findAllSubsets(set, 0);
		ArrayList<ArrayList<Integer>> allSubSets2 = constructSubsets.findAllSubsetsNoRecursion(set);
		for (ArrayList<Integer> s : allSubSets) {
			System.out.println(count + ": " + s);
			count++;
		}
		System.out.println( "--------------------------");
		count = 1;
		for (ArrayList<Integer> s : allSubSets2) {
			System.out.println(count + ": " + s);
			count++;
		}
	}

}
