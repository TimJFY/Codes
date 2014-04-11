package elementary_Data_Structures_Trees.BinaryTrees;

import java.util.ArrayList;

public class DeepestTrails {
	private ArrayList<ArrayList<Node>> deepestTrails = new ArrayList<ArrayList<Node>>();
	private ArrayList<Node> currentTrail = new ArrayList<Node>();

	public void findAllDeepestTrails(Node root) {
		if (root != null) {
			currentTrail.add(root);
			findAllDeepestTrails(root.getLeftChild());
			findAllDeepestTrails(root.getRightChild());
			currentTrail.remove(currentTrail.size() - 1);
		} else {
			// reach the next layer(virtual layer) of a leaf
			if (deepestTrails.size() != 0) {
				if (deepestTrails.get(0).size() < currentTrail.size()) {
					deepestTrails.removeAll(deepestTrails);
					ArrayList<Node> currentTrailCopy = (ArrayList<Node>) currentTrail
							.clone();
					deepestTrails.add(currentTrailCopy);
				} else if (deepestTrails.get(0).size() == currentTrail.size()) {
					// if we find a different trail which has the same length as
					// the longest one, we record it.
					// 'different trails' means the ending nodes are different
					Node currentLast = currentTrail
							.get(currentTrail.size() - 1);
					ArrayList<Node> lastTrail = deepestTrails.get(deepestTrails
							.size() - 1);
					Node previousLast = lastTrail.get(lastTrail.size() - 1);

					if (currentLast != previousLast) {
						ArrayList<Node> currentTrailCopy = (ArrayList<Node>) currentTrail
								.clone();
						deepestTrails.add(currentTrailCopy);
					}
				}
			} else {
				ArrayList<Node> currentTrailCopy = (ArrayList<Node>) currentTrail
						.clone();
				deepestTrails.add(currentTrailCopy);
			}
		}
	}

	public void showAllDeepestTrails() {
		for (int i = 0; i < deepestTrails.size(); i++) {
			System.out.println("Trail: " + (i + 1));
			for (int j = 0; j < deepestTrails.get(i).size(); j++) {
				System.out.print("Node(value = "
						+ deepestTrails.get(i).get(j).getValue() + ") ");
				if (j + 1 < deepestTrails.get(i).size()) {
					if (deepestTrails.get(i).get(j).getLeftChild() == deepestTrails
							.get(i).get(j + 1)) {
						System.out.print("Go_left ");
					} else {
						System.out.print("Go_right ");
					}
				}

			}
			System.out.println();
			System.out.println("------------------------");
		}
	}

	public static void main(String[] args) {
		int[] values = new int[] { 6, 4, 10, -1, -1, 8, -1, -1, 12, 16, -1, -1,
				14, -1, -1 };

		// build the binary tree virtual object
		BinaryTree binaryTree = new BinaryTree();
		// the original binary tree
		Node root = new Node();
		root = binaryTree.crateBinaryTree(values, root, -1);

		DeepestTrails deepestTrails = new DeepestTrails();
		deepestTrails.findAllDeepestTrails(root);
		deepestTrails.showAllDeepestTrails();
	}

}
