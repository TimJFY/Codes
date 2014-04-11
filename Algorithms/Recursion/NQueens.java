package Recursion;

/* 
 * 8.8
 * Write an algorithm to print all ways of arranging eight queens on a chess board 
 * so that none of them share the same row, column or diagonal
 */

public class NQueens {

	// valid column position for each row
	int columnIndexes[];
	int size;
	String[][] chessBorad;
	int solutions = 0;

	public NQueens(int n) {
		this.columnIndexes = new int[n];
		this.size = n;
		this.chessBorad = new String[n][n];
		clearBoard();
	}

	public void placeQueen(int rowIndex) {
		if (rowIndex == size) {
			printBoard();
			solutions++;
			System.out.println("Solution " + solutions + " Completed");
			System.out.println();
			// restore the chess board to initial state
			clearBoard();
			return;
		}
		for (int i = 0; i < size; i++) {
			// as to the current row, try every column
			columnIndexes[rowIndex] = i;
			if (validate(rowIndex)) {
				placeQueen(rowIndex + 1);
			}
		}

	}

	public boolean validate(int rowIndex) {
		for (int j = 0; j < rowIndex; j++) {
			int difference = Math.abs(columnIndexes[j]
					- columnIndexes[rowIndex]);
			if (difference == 0 || difference == rowIndex - j) {
				return false;
			}
		}
		return true;
	}

	public void printBoard() {
		int i = 0;
		for (int j : columnIndexes) {
			chessBorad[i][j] = "*";
			i++;
		}

		for (int row = 0; row < columnIndexes.length; row++) {
			for (int column = 0; column < columnIndexes.length; column++) {
				System.out.print(chessBorad[row][column] + " ");
			}
			System.out.print('\n');
		}
	}

	public void clearBoard() {
		for (int row = 0; row < columnIndexes.length; row++) {
			for (int column = 0; column < columnIndexes.length; column++) {
				chessBorad[row][column] = "-";
			}
		}
	}

	public static void main(String[] args) {
		int n = 8;
		NQueens nQueens = new NQueens(n);
		nQueens.placeQueen(0);
	}
}
