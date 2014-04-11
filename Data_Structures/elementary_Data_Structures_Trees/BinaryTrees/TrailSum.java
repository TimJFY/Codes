package elementary_Data_Structures_Trees.BinaryTrees;

import java.util.ArrayList;

public class TrailSum {

	private ArrayList<ArrayList<Node>> allValidTrails = new ArrayList<ArrayList<Node>>();

	public void testTrail(Node head, int sum, ArrayList<Node> completeTrail,
			int layer) {
		if (head == null) {
			return;
		}

		completeTrail.add(head);
		int temp = sum;
		for (int i = layer; i >= 0; i--) {
			temp -= completeTrail.get(i).getValue();
			// valid
			if (temp == 0) {
				// record all the node along the valid trail
				ArrayList<Node> validTrail = new ArrayList<Node>();
				for (int j = i; j <= layer; j++) {
					validTrail.add(completeTrail.get(j));
				}
				allValidTrails.add(validTrail);
			}
		}

		ArrayList<Node> completeTrailLeft = ((ArrayList<Node>) completeTrail
				.clone());
		ArrayList<Node> completeTrailRight = (ArrayList<Node>) completeTrail
				.clone();

		// go to next layer
		layer++;
		
		testTrail(head.getLeftChild(), sum, completeTrailLeft, layer);
		testTrail(head.getRightChild(), sum, completeTrailRight, layer);
	}

	public void showValidTrails() {
		System.out.println("All Valid Trails:");
		System.out.println("----------------------------");
		for (int i = 0; i < allValidTrails.size(); i++) {
			ArrayList<Node> oneValidTrail = allValidTrails.get(i);
			for (int j = 0; j < oneValidTrail.size(); j++) {
				System.out.print("Node(value = "
						+ oneValidTrail.get(j).getValue() + ") ");
				if (j + 1 < oneValidTrail.size()) {
					if (oneValidTrail.get(j).getLeftChild() == oneValidTrail
							.get(j + 1)) {
						System.out.print("Go_left ");
					} else {
						System.out.print("Go_right ");
					}
				}
			}
			System.out.println();
			System.out.println("----------------------------");
		}
	}

	public static void main(String[] args) {

		int[] values = new int[] { 6, 4, 1, -11, -1000, -1000, -1000, -10,
				-1000, -1000, -6, 0, -1000, -1000, 14, -1000, -1000 };

		// build the binary tree virtual object
		BinaryTree binaryTree = new BinaryTree();
		// the original binary tree
		Node root = new Node();
		root = binaryTree.crateBinaryTree(values, root, -1000);

		TrailSum trailSum = new TrailSum();
		ArrayList<Node> trail = new ArrayList<Node>();
		trailSum.testTrail(root, 0, trail, 0);
		trailSum.showValidTrails();

	}
}
