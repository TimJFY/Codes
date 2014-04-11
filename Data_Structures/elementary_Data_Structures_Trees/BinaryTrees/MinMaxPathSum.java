package elementary_Data_Structures_Trees.BinaryTrees;

import java.util.ArrayList;

public class MinMaxPathSum {

	static int minSum = Integer.MAX_VALUE;
	static int maxSum = Integer.MIN_VALUE;
	static ArrayList<Integer> minSumPath = new ArrayList<Integer>();
	static ArrayList<Integer> maxSumPath = new ArrayList<Integer>();

	public void min_maxPathSum(Node root, int curSum,
			ArrayList<Integer> curMinPath, ArrayList<Integer> curMaxPath) {
		curSum += root.getValue();
		curMinPath.add(root.getValue());
		curMaxPath.add(root.getValue());
		if (root.getLeftChild() == null && root.getRightChild() == null) {
			if (minSum > curSum) {
				minSumPath = curMinPath;
				minSum = curSum;
			}
			if (maxSum < curSum) {
				maxSumPath = curMaxPath;
				maxSum = curSum;
			}
		} else {
			if (root.getLeftChild() != null) {
				min_maxPathSum(root.getLeftChild(), curSum,
						(ArrayList<Integer>) curMinPath.clone(),
						(ArrayList<Integer>) curMaxPath.clone());
			}
			if (root.getRightChild() != null) {
				min_maxPathSum(root.getRightChild(), curSum,
						(ArrayList<Integer>) curMinPath.clone(),
						(ArrayList<Integer>) curMaxPath.clone());
			}
		}
	}

	public static void main(String[] args) {
		int[] values = new int[] { 20, 5, 7, -1000, -1000, -2, -1000, -1000, 3,
				-18, -1, 50, -1000, -1000, 22, -1000, -1000, 19, -1000, -1000,
				-1000 };

//		int[] values = new int[] { 6, -8, 4, -1000, -1000, 5, -1000, -1000, 3,
//				-1000, -1000 };

		// build the binary tree virtual object
		BinaryTree binaryTree = new BinaryTree();
		// the original binary tree
		Node root = new Node();
		root = binaryTree.crateBinaryTree(values, root, -1000);

		MinMaxPathSum pathSum = new MinMaxPathSum();
		pathSum.min_maxPathSum(root, 0, new ArrayList<Integer>(),
				new ArrayList<Integer>());
		System.out.println("MinPathSum = " + minSum);
		System.out.println("MinPath: " + minSumPath);

		System.out.println("MaxPathSum = " + maxSum);
		System.out.println("MaxPath: " + maxSumPath);

	}
}
