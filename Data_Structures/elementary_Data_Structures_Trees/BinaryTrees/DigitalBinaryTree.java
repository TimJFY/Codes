package elementary_Data_Structures_Trees.BinaryTrees;

import java.util.ArrayList;
import java.util.Random;

public class DigitalBinaryTree {

	// Data structure used to store all the node elements of the binaryTree
	private ArrayList<Node> BinaryTree = new ArrayList<Node>();

	/**
	 * Construct the binary tree, (n-->2n, n-->2n+1)
	 */
	public DigitalBinaryTree(int nodeAmount) {
		if (nodeAmount <= 0) {
			System.out
					.println("False, can only creat with positive number of nodes");
			return;
		}
		if (nodeAmount == 1) {
			Node oneNode = new Node(1);
			BinaryTree.add(oneNode);
		} else {
			for (int i = 1; i <= nodeAmount; i++) {
				Node oneNode = new Node(i);
				BinaryTree.add(oneNode);
			}
			int j = nodeAmount / 2;
			// Generate the cost of each path Randomly
			Random randomCost = new Random();
			for (int k = 1; k < j; k++) {
				// For each parent node has a value i, the values of its
				// left-child node and right-child node should be (2*i) and
				// (2*+1) respectively
				BinaryTree.get(k - 1).setLeftChild(BinaryTree.get(2 * k - 1));
				// The cost of each path will be a random integer value between
				// [0,15]
				BinaryTree.get(k - 1)
						.setToLeftChildCost(randomCost.nextInt(16));
				BinaryTree
						.get(k - 1)
						.getLeftChild()
						.setTotalReachCost(
								BinaryTree.get(k - 1).getToLeftChildCost()
										+ BinaryTree.get(k - 1)
												.getTotalReachCost());
				BinaryTree.get(k - 1).setRightChild(BinaryTree.get(2 * k));
				BinaryTree.get(k - 1).setToRightChildCost(
						randomCost.nextInt(16));
				BinaryTree
						.get(k - 1)
						.getRightChild()
						.setTotalReachCost(
								BinaryTree.get(k - 1).getToRightChildCost()
										+ BinaryTree.get(k - 1)
												.getTotalReachCost());
			}
			// If the total number of nodes in the tree is even, then the last
			// node that has child node must only have a left-child node
			if (nodeAmount % 2 == 0) {
				BinaryTree.get(j - 1).setLeftChild(
						BinaryTree.get(nodeAmount - 1));
				BinaryTree.get(j - 1)
						.setToLeftChildCost(randomCost.nextInt(16));
				BinaryTree
						.get(j - 1)
						.getLeftChild()
						.setTotalReachCost(
								BinaryTree.get(j - 1).getToLeftChildCost()
										+ BinaryTree.get(j - 1)
												.getTotalReachCost());
			} else {
				BinaryTree.get(j - 1).setLeftChild(
						BinaryTree.get(nodeAmount - 2));
				BinaryTree.get(j - 1)
						.setToLeftChildCost(randomCost.nextInt(16));
				BinaryTree
						.get(j - 1)
						.getLeftChild()
						.setTotalReachCost(
								BinaryTree.get(j - 1).getToLeftChildCost()
										+ BinaryTree.get(j - 1)
												.getTotalReachCost());
				BinaryTree.get(j - 1).setRightChild(
						BinaryTree.get(nodeAmount - 1));
				BinaryTree.get(j - 1).setToRightChildCost(
						randomCost.nextInt(16));
				BinaryTree
						.get(j - 1)
						.getRightChild()
						.setTotalReachCost(
								BinaryTree.get(j - 1).getToRightChildCost()
										+ BinaryTree.get(j - 1)
												.getTotalReachCost());
			}
		}
	}

	public DigitalBinaryTree() {
	}

	/**
	 * Compute the depth of the tree
	 */
	public int treeDepth() {
		int nodeAmount = this.BinaryTree.size();
		int depth = (int) (Math.log(nodeAmount) / Math.log(2)) + 1;
		return depth;
	}

	/**
	 * Trim a certain binary tree according to the depth needed
	 */
	public DigitalBinaryTree treeTrim(DigitalBinaryTree binaryTree, int depth) {
		int index = (int) Math.pow(2, (depth - 1));
		int nodesAmount = 0;
		if ((int) Math.pow(2, depth) < binaryTree.getBinaryTree().size()) {
			nodesAmount = (int) Math.pow(2, depth) - 1;
		} else {
			nodesAmount = binaryTree.getBinaryTree().size();
		}
		DigitalBinaryTree newBinaryTree = new DigitalBinaryTree();
		ArrayList<Node> newBinaryTrees = newBinaryTree.getBinaryTree();
		// Construct a new tree whose attributes(except for the totalReachCost,
		// useless ) are exactly identical with the former tree
		if (nodesAmount == 1) {
			Node oneNode = new Node(1);
			newBinaryTrees.add(oneNode);
		} else {
			for (int i = 1; i <= nodesAmount; i++) {
				Node oneNode = new Node(i);
				newBinaryTrees.add(oneNode);
			}
			int j = nodesAmount / 2;

			for (int k = 1; k < j; k++) {
				newBinaryTrees.get(k - 1).setLeftChild(
						newBinaryTrees.get(2 * k - 1));
				newBinaryTrees.get(k - 1).setToLeftChildCost(
						binaryTree.getBinaryTree().get(k - 1)
								.getToLeftChildCost());
				newBinaryTrees.get(k - 1).setRightChild(
						newBinaryTrees.get(2 * k));
				newBinaryTrees.get(k - 1).setToRightChildCost(
						binaryTree.getBinaryTree().get(k - 1)
								.getToRightChildCost());
			}
			if (nodesAmount % 2 == 0) {
				newBinaryTrees.get(j - 1).setLeftChild(
						newBinaryTrees.get(nodesAmount - 1));
				newBinaryTrees.get(j - 1).setToLeftChildCost(
						binaryTree.getBinaryTree().get(j - 1)
								.getToLeftChildCost());
			} else {
				newBinaryTrees.get(j - 1).setLeftChild(
						newBinaryTrees.get(nodesAmount - 2));
				newBinaryTrees.get(j - 1).setToLeftChildCost(
						binaryTree.getBinaryTree().get(j - 1)
								.getToLeftChildCost());
				newBinaryTrees.get(j - 1).setRightChild(
						newBinaryTrees.get(nodesAmount - 1));
				newBinaryTrees.get(j - 1).setToRightChildCost(
						binaryTree.getBinaryTree().get(j - 1)
								.getToRightChildCost());
			}
		}

		// Trim all nodes that beyond the depth needed
		while (index < (int) Math.pow(2, depth)
				&& index < binaryTree.getBinaryTree().size()) {
			newBinaryTree.getBinaryTree().get(index - 1).setLeftChild(null);
			newBinaryTree.getBinaryTree().get(index - 1).setRightChild(null);
			index++;
		}
		return newBinaryTree;
	}

	/**
	 * Output the binary tree structure
	 */
	public void display() {
		System.out.println("Here is the binary tree:");
		for (int i = 0; i < BinaryTree.size(); i++) {
			System.out.println(BinaryTree.get(i).toString());
		}
		System.out.println();
	}

	public ArrayList<Node> getBinaryTree() {
		return BinaryTree;
	}

	public void setBinaryTree(ArrayList<Node> binaryTree) {
		BinaryTree = binaryTree;
	}
}
