package elementary_Data_Structures_Stacks_And_Queues;

import java.util.Stack;

public class HanoiTowers {

	Tower[] towers;

	public HanoiTowers(int n) {

		Tower tower1 = new Tower(1);
		// put n plates into the first tower
		for (int i = n - 1; i >= 0; i--) {
			tower1.addDisk(i);
		}

		Tower tower2 = new Tower(2);
		Tower tower3 = new Tower(3);

		towers = new Tower[3];
		towers[0] = tower1;
		towers[1] = tower2;
		towers[2] = tower3;
	}

	public static void main(String args[]) {
		int n = 5;

		HanoiTowers hanoiTowers = new HanoiTowers(n);
		System.out.print("Inintial State: Tower 1 :");
		for (int j = 0; j <= hanoiTowers.towers[0].getDisks().size() - 1; j++) {
			System.out.print(" " + hanoiTowers.towers[0].getDisks().get(j));
		}
		System.out.println();
		System.out.println("------------------");
		hanoiTowers.towers[0].moveDisks(n, hanoiTowers.towers[2],
				hanoiTowers.towers[1], hanoiTowers.towers);
	}

	class Tower {
		int index;
		Stack<Integer> disks;

		public Tower(int index) {
			this.index = index;
			this.disks = new Stack<Integer>();
		}

		public int getIndex() {
			return index;
		}

		public Stack<Integer> getDisks() {
			return disks;
		}

		public void moveDisks(int num, Tower destination, Tower buffer,
				Tower[] allTowers) {
			if (num > 0) {
				// move first (n-1) plates from origin to buffer tower
				// using destination tower as 'buffer'
				moveDisks(num - 1, buffer, destination, allTowers);
				// move the last(largest) plate from origin to destination tower
				moveTop(destination, allTowers);
				// move first (n-1) plates from buffer to destination tower
				// using origin as 'buffer'
				buffer.moveDisks(num - 1, destination, this, allTowers);

			}
		}

		public void moveTop(Tower destination, Tower[] allTowers) {
			int top = this.getDisks().pop();
			destination.addDisk(top);
			System.out.println("Move disk " + top + " form Tower"
					+ this.getIndex() + " to Tower" + destination.getIndex());
			showInfo(allTowers);
			System.out.println("---------------");
		}

		public void addDisk(int size) {
			if (!this.getDisks().isEmpty() && this.getDisks().peek() < size) {
				System.out.println("Error when placing disk : " + size);
			} else {
				this.getDisks().push(size);
			}

		}

		public void showInfo(Tower[] allTowers) {
			for (int i = 0; i < allTowers.length; i++) {
				System.out.print("Contents in tower " + allTowers[i].getIndex()
						+ " :");
				for (int j = 0; j <= allTowers[i].getDisks().size() - 1; j++) {
					System.out.print(" " + allTowers[i].getDisks().get(j));
				}
				System.out.println();
			}

		}
	}
}
