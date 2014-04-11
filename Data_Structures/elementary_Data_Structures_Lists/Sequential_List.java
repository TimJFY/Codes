package elementary_Data_Structures_Lists;

import java.util.ArrayList;

public class Sequential_List {
	private int list_size;
	private ArrayList<Integer> seq_List;

	public Sequential_List() {
		list_size = 0;
		seq_List = new ArrayList<Integer>();
	}

	public void insertElement(int index, int value) {
		int offset;
		if (index > list_size || index < 0) {
			System.out.println("Invalid Parameter");
		} else {
			for (offset = list_size - 1; offset >= index; offset--) {
				seq_List.set(offset + 1, seq_List.get(offset));
			}
			seq_List.set(index, value);
			list_size++;
		}
	}

	public void removeElement(int index) {
		int offset;
		if (index >= list_size || index < 0) {
			System.out.println("Invalid Parameter");
		} else {
			for (offset = index + 1; offset < list_size; offset++) {
				seq_List.set(offset - 1, seq_List.get(offset));
			}

			list_size--;
		}
	}

	public void setElement(int index, int value) {
		this.seq_List.set(index, value);
	}

	public int getElement(int index) {
		return this.seq_List.get(index);
	}

}
