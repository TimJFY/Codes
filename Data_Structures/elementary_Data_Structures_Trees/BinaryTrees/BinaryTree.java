package elementary_Data_Structures_Trees.BinaryTrees;

import java.util.LinkedList;
import java.util.Queue;
import java.util.Stack;

public class BinaryTree {

	public static int startIndex = 0;

	public Node crateBinaryTree(int[] values, Node root, int nullNum) {
		int value = values[startIndex++];
		// nullNum means no child
		if (value == nullNum) {
			root = null;
		} else {
			root = new Node(value, null, null);
			root.setLeftChild(crateBinaryTree(values, root.getLeftChild(),
					nullNum));
			root.setRightChild(crateBinaryTree(values, root.getRightChild(),
					nullNum));
		}
		return root;
	}

	public void preOrderTraverse(Node root) {
		if (!(root == null)) {
			System.out.print(" " + root.getValue() + " ");
			preOrderTraverse(root.getLeftChild());
			preOrderTraverse(root.getRightChild());
		}
	}

	public void inOrderTraverse(Node root) {
		if (!(root == null)) {
			inOrderTraverse(root.getLeftChild());
			System.out.print(" " + root.getValue() + " ");
			inOrderTraverse(root.getRightChild());
		}
	}

	public void postOrderTraverse(Node root) {
		if (!(root == null)) {
			postOrderTraverse(root.getLeftChild());
			postOrderTraverse(root.getRightChild());
			System.out.print(" " + root.getValue() + " ");
		}
	}

	public void breadthFirstSearch(Node root) {
		if (!(root == null)) {
			Queue<Node> bts_Queue = new LinkedList<Node>();
			bts_Queue.add(root);
			while (!bts_Queue.isEmpty()) {
				Node currentNode = bts_Queue.poll();
				System.out.print(" " + currentNode.getValue() + " ");
				if (currentNode.getLeftChild() != null) {
					bts_Queue.offer(currentNode.getLeftChild());
				}
				if (currentNode.getRightChild() != null) {
					bts_Queue.offer(currentNode.getRightChild());
				}
			}
		} else
			System.out.print("Empty Tree");
	}

	public void depthFirstSearch(Node root) {
		if (!(root == null)) {
			Stack<Node> bts_Queue = new Stack<Node>();
			bts_Queue.add(root);
			while (!bts_Queue.isEmpty()) {
				Node currentNode = bts_Queue.pop();
				System.out.print(" " + currentNode.getValue() + " ");
				if (currentNode.getRightChild() != null) {
					bts_Queue.push(currentNode.getRightChild());
				}
				// left child will be read first
				if (currentNode.getLeftChild() != null) {
					bts_Queue.push(currentNode.getLeftChild());
				}
			}
		} else
			System.out.print("Empty Tree");
	}

	public void refresh() {
		startIndex = 0;
	}

}
