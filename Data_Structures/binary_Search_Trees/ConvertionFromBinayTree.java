package binary_Search_Trees;

import java.util.LinkedList;
import java.util.Queue;

import elementary_Data_Structures_Trees.BinaryTrees.BinaryTree;
import elementary_Data_Structures_Trees.BinaryTrees.Node;

public class ConvertionFromBinayTree {

	public Node binaryTreeToBST(Node root) {
		// create a new tree, need additional space
		Node searchTreeRoot = null;
		if (root != null) {
			// use BFS to insert all nodes into a new BST
			Queue<Node> bts_Queue = new LinkedList<Node>();
			bts_Queue.add(root);
			while (!bts_Queue.isEmpty()) {
				Node currentNode = bts_Queue.poll();
				searchTreeRoot = insertNode(searchTreeRoot,
						currentNode.getValue());
				if (currentNode.getLeftChild() != null) {
					bts_Queue.offer(currentNode.getLeftChild());
				}
				if (currentNode.getRightChild() != null) {
					bts_Queue.offer(currentNode.getRightChild());
				}
			}

		}
		return searchTreeRoot;
	}

	public Node binaryTreeToBSTinPlace(Node root, Node bstRoot) {
		// do not use additional space
		if (!(root == null)) {
			// use post-order-traverse to link all nodes to the root
			bstRoot = binaryTreeToBSTinPlace(root.getLeftChild(), bstRoot);
			bstRoot = binaryTreeToBSTinPlace(root.getRightChild(), bstRoot);
			// link an old to the root according to the rules of BST
			bstRoot = reInsertNode(bstRoot, root);
		}
		return bstRoot;
	}

	// insert node to BST , inplace
	public Node reInsertNode(Node bstRoot, Node leaf) {
		if (bstRoot != null) {
			if (leaf.getValue() < bstRoot.getValue()) {
				bstRoot.setLeftChild(reInsertNode(bstRoot.getLeftChild(), leaf));
			} else if (leaf.getValue() > bstRoot.getValue()) {
				bstRoot.setRightChild(reInsertNode(bstRoot.getRightChild(),
						leaf));
			}
		} else {
			leaf.setLeftChild(null);
			leaf.setRightChild(null);
			// link to previous node
			bstRoot = leaf;
		}
		return bstRoot;
	}

	// insert a node to BST, create a new tree
	public Node insertNode(Node leaf, int nodeValue) {
		if (leaf == null) {
			// new a node
			leaf = new Node(nodeValue, null, null);
		} else if (nodeValue < leaf.getValue()) {
			leaf.setLeftChild(insertNode(leaf.getLeftChild(), nodeValue));
		} else if (nodeValue > leaf.getValue()) {
			leaf.setRightChild(insertNode(leaf.getRightChild(), nodeValue));
		}
		return leaf;
	}

	public static void main(String[] args) {

		// serialization of a binary tree
		int[] values = new int[] { 6, 4, 10, 5, -1, -1, -1, 8, -1, -1, 12, 16,
				-1, -1, 14, -1, -1 };

		// build the binary tree virtual object
		BinaryTree binaryTree = new BinaryTree();
		// root of the original binary tree
		Node root = new Node();
		// construct the tree by its pre-order traversal
		root = binaryTree.crateBinaryTree(values, root, -1);
		System.out.println("Original Binary Tree Pre-Order-Traversal:");
		binaryTree.preOrderTraverse(root);
		System.out.println();
		System.out.println("--------------------------------------------");
		// binaryTree.inOrderTraverse(root);
		// binaryTree.postOrderTraverse(root);
		// System.out.println("Breadth First Traversal:");
		// binaryTree.breadthFirstSearch(root);
		// System.out.println();
		// System.out.println("Depth First Traversal:");
		// binaryTree.depthFirstSearch(root);
		// System.out.println();

		ConvertionFromBinayTree binayTreeConvert = new ConvertionFromBinayTree();
		Node searchTreeRoot = binayTreeConvert.binaryTreeToBST(root);

		System.out.println("BST, a new tree");
		System.out.println("BST Pre-Order-Traversal:");
		binaryTree.preOrderTraverse(searchTreeRoot);
		System.out.println();
		System.out.println("--------------------------------------------");

		// only allocate space for the new root(value not changed)
		Node bstRoot = new Node(root.getValue());
		bstRoot = binayTreeConvert.binaryTreeToBSTinPlace(root, bstRoot);

		System.out.println("Another BST, inplace converted");
		System.out.println("Pre-Order-Traversal:");
		binaryTree.preOrderTraverse(bstRoot);
		System.out.println();
		System.out.println("Post-Order-Traversal:");
		binaryTree.postOrderTraverse(bstRoot);
		System.out.println();
		System.out.println("BFS Traversal:");
		binaryTree.breadthFirstSearch(bstRoot);
		System.out.println();
		System.out.println("--------------------------------------------");

		System.out
				.println("Original Binary Tree Pre-Order-Traversal: (destroyed by inplace convertion)");
		binaryTree.preOrderTraverse(root);
	}
}
