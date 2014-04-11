package elementary_Data_Structures_Lists;

public class Linked_List {

	public Node createLinkedList(int[] valueArray) {
		if (valueArray.length == 0) {
			return null;
		}
		Node head = new Node(valueArray[0], null);
		Node pointer = head;
		for (int i = 1; i < valueArray.length; i++) {
			Node element = new Node(valueArray[i], null);
			pointer.setNext(element);
			pointer = pointer.getNext();
		}

		return head;
	}

	public void show(Node head) {
		while (head != null) {
			System.out.print(head.getValue());
			head = head.getNext();
			if (head != null) {
				System.out.print(" ---> ");
			}
		}
		System.out.println();
	}

}
