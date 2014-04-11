package elementary_Data_Structures_Trees.BinaryTrees;

import java.util.ArrayList;

public class LowsetComAncestor {

	ArrayList<Node> pathA = new ArrayList<Node>();
	ArrayList<Node> pathA_ = new ArrayList<Node>();
	ArrayList<Node> pathB = new ArrayList<Node>();
	ArrayList<Node> pathB_ = new ArrayList<Node>();

	// find the path from the root to a target node
	public boolean findPath(Node root, Node c, int index) {
		Node current = root;
		boolean isAns = false;
		if (current.getValue() == c.getValue()) {
			isAns = true;
		} else {
			if (current.getLeftChild() != null && !isAns) {
				isAns = findPath(current.getLeftChild(), c, index);
			}
			if (current.getRightChild() != null && !isAns) {
				isAns = findPath(current.getRightChild(), c, index);
			}
		}
		// construct the path in a bottom-up manner
		if (isAns && index == 1)
			pathA.add(current);
		if (isAns && index == 2)
			pathB.add(current);
		return isAns;
	}

	// find the path from the root to a target node
	public boolean findPath2(Node root, Node c, int index,
			ArrayList<Node> curPath) {
		// construct the path in a top-down manner
		curPath.add(root);
		boolean found = false;
		if (root.getValue() == c.getValue()) {
			if (index == 1) {
				pathA_ = curPath;
			} else {
				pathB_ = curPath;
			}
			found = true;
		} else {
			if (root.getLeftChild() != null && !found) {
				found = findPath2(root.getLeftChild(), c, index,
						(ArrayList<Node>) curPath.clone());
			}
			if (root.getRightChild() != null && !found) {
				found = findPath2(root.getRightChild(), c, index,
						(ArrayList<Node>) curPath.clone());
			}
		}
		return found;
	}

	public void showPaths() {
		System.out.println("Bottom-Up");
		for (int i = 0; i < pathA.size(); i++) {
			System.out.print(pathA.get(i).getValue());
			if (i + 1 < pathA.size()) {
				System.out.print(" --> ");
			} else {
				System.out.print(" : 1st node");
			}
		}
		System.out.println();
		for (int i = 0; i < pathB.size(); i++) {
			System.out.print(pathB.get(i).getValue());
			if (i + 1 < pathB.size()) {
				System.out.print(" --> ");
			} else {
				System.out.print(" : 2nd node");
			}
		}
		System.out.println();
		System.out.println("Top-Down");
		for (int i = 0; i < pathA_.size(); i++) {
			System.out.print(pathA_.get(i).getValue());
			if (i + 1 < pathA_.size()) {
				System.out.print(" --> ");
			} else {
				System.out.print(" : 1st node");
			}
		}
		System.out.println();
		for (int i = 0; i < pathB_.size(); i++) {
			System.out.print(pathB_.get(i).getValue());
			if (i + 1 < pathB_.size()) {
				System.out.print(" --> ");
			} else {
				System.out.print(" : 2nd node");
			}
		}
		System.out.println();

	}

	public void findLowestCommonAncestor() {
		if (pathA.size() == 0 || pathA_.size() == 0) {
			System.out.println("1st node does not exist in the tree");
			return;
		} else if (pathB.size() == 0 || pathB_.size() == 0) {
			System.out.println("2nd node does not exist in the tree");
			return;
		}
		int minLength = Math.min(pathA_.size(), pathB_.size()) - 1;
		int commonNode = 0;
		while (commonNode <= minLength) {
			
			if (pathA_.get(commonNode) != pathB_.get(commonNode)) {
				System.out.println("The lowset common ancestor is: "
						+ pathA_.get(commonNode - 1));
				break;
			}
			if (commonNode == minLength) {
				if (pathA_.size() == pathB_.size()) {
					System.out
							.println("Same node, The lowset common ancestor is: "
									+ pathA_.get(commonNode - 1));
				} else {
					System.out
							.println("The lowset common ancestor is: "
									+ pathA_.get(commonNode));
				}
			}
			
			commonNode++;
		}
	}

	public static void main(String[] args) {

		Node n31 = new Node(8, null, null);
		Node n32;
		Node n33 = new Node(10, null, null);
		Node n34;
		Node n35;
		Node n36;
		Node n37 = new Node(14, null, null);
		Node n38 = new Node(15, null, null);

		Node n21 = new Node(4, n31, null);
		Node n22 = new Node(5, n33, null);
		Node n23 = new Node(6, null, null);
		Node n24 = new Node(7, n37, n38);

		Node n11 = new Node(2, n21, n22);
		Node n12 = new Node(3, n23, n24);

		Node root = new Node(1, n11, n12);

		Node notExist = new Node(100, null, null);

		LowsetComAncestor lowComAns = new LowsetComAncestor();
		lowComAns.findPath(root, n21, 1);
		lowComAns.findPath(root, n33, 2);
		lowComAns.findPath2(root, n21, 1, new ArrayList<Node>());
		lowComAns.findPath2(root, n33, 2, new ArrayList<Node>());
		lowComAns.showPaths();
		lowComAns.findLowestCommonAncestor();
	}

}