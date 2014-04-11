package elementary_Data_Structures_Lists;

public class Node {
	private int value;
	Node next;

	public Node(int value, Node next) {
		this.value = value;
		this.next = next;
	}

	public Node(int value) {
		this.value = value;
		this.next = null;
	}

	public int getValue() {
		return value;
	}

	public void setValue(int value) {
		this.value = value;
	}

	public Node getNext() {
		return next;
	}

	public void setNext(Node next) {
		this.next = next;
	}

}
