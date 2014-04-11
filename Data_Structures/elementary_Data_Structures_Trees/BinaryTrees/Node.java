package elementary_Data_Structures_Trees.BinaryTrees;

public class Node {
	// Attribute NodeValue
	private int value = 0;
	// Attribute leftChild
	private Node leftChild;
	// Attribute rightChild
	private Node rightChild;
	// Attribute the cost from parent Node to its LeftChild Node
	private int toLeftChildCost = 0;
	// Attribute the cost from parent Node to its RightChild Node
	private int toRightChildCost = 0;
	// Attribute the cost from root to the current Node
	private int totalReachCost = 0;

	public Node(int value) {
		this.value = value;
		this.leftChild = null;
		this.rightChild = null;
	}

	public Node(int value, Node leftChild, Node rightChild) {
		this.value = value;
		this.leftChild = leftChild;
		this.rightChild = rightChild;
	}

	public Node() {
	}

	public int getValue() {
		return value;
	}

	public void setValue(int value) {
		this.value = value;
	}

	public Node getLeftChild() {
		return leftChild;
	}

	public void setLeftChild(Node leftChild) {
		this.leftChild = leftChild;
	}

	public Node getRightChild() {
		return rightChild;
	}

	public void setRightChild(Node rightChild) {
		this.rightChild = rightChild;
	}

	public int getToLeftChildCost() {
		return toLeftChildCost;
	}

	public void setToLeftChildCost(int toLeftChildCost) {
		this.toLeftChildCost = toLeftChildCost;
	}

	public int getToRightChildCost() {
		return toRightChildCost;
	}

	public void setToRightChildCost(int toRightChildCost) {
		this.toRightChildCost = toRightChildCost;
	}

	public int getTotalReachCost() {
		return totalReachCost;
	}

	public void setTotalReachCost(int totalReachCost) {
		this.totalReachCost = totalReachCost;
	}

	/**
	 * Output the node information(All Attributes)
	 */
	public String toString() {
		String s = "parentNode value: " + this.value + " totalReachCost= "
				+ this.totalReachCost;
		if (null == this.getLeftChild()) {
			s += " ---->No Child ";
		} else if (null == this.getRightChild()) {
			s += " ---->Only leftChild " + this.leftChild.getValue()
					+ " Cost = " + this.getToLeftChildCost();
		} else {
			s += " ---->leftChild value: " + this.leftChild.getValue()
					+ " Cost = " + this.getToLeftChildCost()
					+ " ---->rightChild value: " + this.rightChild.getValue()
					+ " Cost = " + this.getToRightChildCost();
		}
		return s;
	}

	/**
	 * Output the tree structure(Relationship among nodes)
	 */
	public String showStructure() {
		if (this.getLeftChild() == null && this.getRightChild() == null) {
			return "(" + this.getValue() + ")";
		} else {
			return "(" + this.getValue() + ", " + this.getLeftChild().showStructure()
					+ ", " + this.getRightChild().showStructure() + ")";
		}
	}
}
