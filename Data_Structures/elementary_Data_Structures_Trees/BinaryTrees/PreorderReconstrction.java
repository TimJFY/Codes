package elementary_Data_Structures_Trees.BinaryTrees;

import java.util.*;

// reconstruct a full binary tree using its Preorder
class PreorderReconstrction {

	private List<Integer> preorder;

	public PreorderReconstrction(String preorderList) {
		preorder = new ArrayList<Integer>();
		for (int i = 0; i < preorderList.length(); i++) {
			preorder.add(new Integer(preorderList.charAt(i) - '0'));
		}
	}

	public Node buildTree() {
		if (preorder.size() > 0) {
			int token = preorder.remove(0);
			if (token == 1) {
				return new Node(token, null, null);
			} else {
				return new Node(token, buildTree(), buildTree());

			}
		} else
			return null;
	}

	public static void main(String args[]) {
		// '0': nodes whose degree are 2, '1': nodes whose degree are 0.
		// necessary conditions for valid input : #'1's = #'0's + 1
		String input = "001101000111011";
		// invalid input, discard some elements
		String input1 = "001110";
		PreorderReconstrction preorder = new PreorderReconstrction(input);
		Node tree = preorder.buildTree();
		System.out.println("Valid:");
		System.out.println(tree.showStructure());
		
		PreorderReconstrction preorder1 = new PreorderReconstrction(input1);
		Node tree1 = preorder1.buildTree();
		System.out.println("The sequence is Invalid:");
		System.out.println(tree1.showStructure());
	}
}
