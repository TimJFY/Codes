package elementary_Data_Structures_Lists;

public class SwapNodesPairs {
	public Node swapPairs(Node head) {
		Node prev = new Node(0);
		Node newhead = prev;
		Node curr = head;
		prev.next = curr;
		while (curr != null && curr.next != null) {
			Node next = curr.next.next;
			prev.next = curr.next;
			curr.next.next = curr;
			curr.next = next;
			prev = curr;
			curr = next;
		}
		return newhead.next;
	}

	public static void main(String[] args) {
		Node d = new Node(1);
		Node c = new Node(2, d);
		Node b = new Node(3, c);
		Node a = new Node(4, b);
		SwapNodesPairs swapNodesPairs = new SwapNodesPairs();
		swapNodesPairs.swapPairs(a);

	}
}
