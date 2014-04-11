package binary_Search_Trees;

import elementary_Data_Structures_Trees.BinaryTrees.BinaryTree;
import elementary_Data_Structures_Trees.BinaryTrees.Node;

public class BSTtoDoubleLinkedList {

	// recursively construct the left subtree and the right subtree
	public Node binarySearchTreeToDoubleLinkedList1(Node root, boolean asRight) {
		// root of the sub tree
		if (root == null) {
			return null;
		}

		Node leftSide = null;
		Node rightSide = null;
		Node linkedListHead = root;

		if (root.getLeftChild() != null) {
			leftSide = binarySearchTreeToDoubleLinkedList1(root.getLeftChild(),
					false);
		}
		// connect the largest node in the left side to the current node
		if (leftSide != null) {
			leftSide.setRightChild(root);
			root.setLeftChild(leftSide);
		}
		if (root.getRightChild() != null) {
			rightSide = binarySearchTreeToDoubleLinkedList1(
					root.getRightChild(), true);
		}
		// connect the smallest node in the right side to the current node
		if (rightSide != null) {
			rightSide.setLeftChild(root);
			root.setRightChild(rightSide);
		}

		// if the current node is the left child of its parent, then the linked
		// list head should be the largest child of the current node
		if (!asRight) {
			while (linkedListHead.getRightChild() != null) {
				linkedListHead = linkedListHead.getRightChild();
			}
		}
		// if the current node is the right child of its parent, then the linked
		// list head should be the smallest child of the current node
		else {
			while (linkedListHead.getLeftChild() != null) {
				linkedListHead = linkedListHead.getLeftChild();
			}
		}
		return linkedListHead;
	}

	// use In-Order-Traverse to construct the linked list
	public Node binarySearchTreeToDoubleLinkedList2(Node root,
			Node linkedListTail) {
		if (root == null) {
			return null;
		}
		if (root.getLeftChild() != null) {
			linkedListTail = binarySearchTreeToDoubleLinkedList2(
					root.getLeftChild(), linkedListTail);
		}
		// put the current node into the double linked list
		root.setLeftChild(linkedListTail);
		if (linkedListTail != null) {
			linkedListTail.setRightChild(root);
		}
		// current node becomes the tail(the last element in the linked list)
		linkedListTail = root;

		if (root.getRightChild() != null) {
			linkedListTail = binarySearchTreeToDoubleLinkedList2(
					root.getRightChild(), linkedListTail);
		}
		return linkedListTail;
	}

	public void showLinkedList(Node head) {
		if (head != null) {
			System.out.println("Double linked List:");
			System.out.println("Forward Travesal:");
			System.out.print(head.getValue());
			while (head.getRightChild() != null) {
				head = head.getRightChild();
				System.out.print(" ---> " + head.getValue());
			}
			System.out.println();
			System.out.println("Backward Travesal:");
			System.out.print(head.getValue());
			while (head.getLeftChild() != null) {
				head = head.getLeftChild();
				System.out.print(" ---> " + head.getValue());
			}
			System.out.println();
		}
	}

	public static void main(String[] args) {
		int[] values = new int[] { 6, 5, 4, -1, -1, -1, 10, 8, -1, -1, 16, 14,
				-1, -1, 12, -1, -1 };
		// build the binary tree virtual object
		BinaryTree binaryTree = new BinaryTree();
		// the original binary search tree
		Node root1 = new Node();
		Node root2 = new Node();
		root1 = binaryTree.crateBinaryTree(values, root1, -1);
		// restore the static index to 0 before create a new tree
		binaryTree.refresh();
		root2 = binaryTree.crateBinaryTree(values, root2, -1);
		System.out.println("Original Binary Search Tree Pre-Order-Traversal:");
		binaryTree.preOrderTraverse(root1);
		System.out.println();
		System.out.println("--------------------------------------------");

		Node linkedListHead1;
		Node linkedListHead2;
		BSTtoDoubleLinkedList bSTtoDoubleLinkedList = new BSTtoDoubleLinkedList();
		// "true" means we will return the list head
		linkedListHead1 = bSTtoDoubleLinkedList
				.binarySearchTreeToDoubleLinkedList1(root1, true);
		bSTtoDoubleLinkedList.showLinkedList(linkedListHead1);
		System.out.println("--------------------------------------------");
		linkedListHead2 = bSTtoDoubleLinkedList
				.binarySearchTreeToDoubleLinkedList2(root2, null);
		bSTtoDoubleLinkedList.showLinkedList(linkedListHead2);
	}

}
