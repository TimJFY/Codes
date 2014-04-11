package elementary_Data_Structures_Lists;

import java.util.Hashtable;

public class RemoveRepet {
	// use hash table, additional space needed
	public Node removeRepeatElements(Node head) {
		Hashtable<Integer, Boolean> distinctElements = new Hashtable<Integer, Boolean>();
		Node pointer = null;
		Node listhHead = head;

		while (head != null) {
			if (distinctElements.containsKey(head.getValue())) {
				pointer.setNext(head.getNext());
			} else {
				distinctElements.put(head.getValue(), true);
				pointer = head;
			}
			head = head.getNext();
		}
		return listhHead;
	}

	// use two indexes, change in place
	public Node removeRepeatElementsInPalce(Node head) {

		Node previous = head;
		if (previous == null) {
			return null;
		}

		Node current = previous.getNext();
		while (current != null) {
			// start from the head each time
			Node pointer = head;
			// we hope there are no duplicates
			boolean noDuplicated = true;
			while (true) {
				if (pointer.getValue() == current.getValue()) {
					noDuplicated = false;
					break;
				}
				if (pointer == previous) {
					break;
				}
				pointer = pointer.getNext();
			}
			// no duplicated element is found(corresponding to current node) 
			// in previous distinct list
			if (noDuplicated) {
				previous = current;
			} else { // drop the current node
				previous.setNext(current.getNext());
			}
			current = current.getNext();
		}
		return head;
	}

	public static void main(String[] args) {
		int[] valueArray = new int[] { 1, 2, 3, 1, 2, 3, -1, 0, 4, 4, 5, 6, 6,
				5, 7, 8, 9 };
		int[] valueArray2 = new int[] { -1, -1, -1, 5, -1, -1, 6, 6, 0, -1, -1 };
		Linked_List linked_List = new Linked_List();
		Node listHead = linked_List.createLinkedList(valueArray);
		Node listHead2 = linked_List.createLinkedList(valueArray2);

		RemoveRepet removeRepet = new RemoveRepet();
		listHead = removeRepet.removeRepeatElements(listHead);
		linked_List.show(listHead);

		listHead2 = removeRepet.removeRepeatElementsInPalce(listHead2);
		linked_List.show(listHead2);

	}
}