package elementary_Data_Structures_Arrays_And_Strings;

/*
 * If an N X N matrix is given, print it in spiral order. 
 */
public class MatrixSpiralRead {

	public void spiralRead(char[][] matrix) {
		int layer = matrix.length / 2;
		// print the whole matrix from outer-most layer to inner-most layer
		for (int i = 0; i <= layer; i++) {
			int start = i;
			int end = matrix.length - 1 - i;
			// if matrix length is an odd number, the most inner layer will be a
			// single grid
			if (start == end) {
				System.out.print(matrix[start][end]);
			} else {
				// left - right
				for (int column = start; column < end; column++) {
					System.out.print(matrix[start][column] + " ");
				}
				// up - down
				for (int row = start; row < end; row++) {
					System.out.print(matrix[row][end] + " ");
				}
				// right - left
				for (int column = end; column > start; column--) {
					System.out.print(matrix[end][column] + " ");
				}
				// down - up
				for (int row = end; row > start; row--) {
					System.out.print(matrix[row][start] + " ");
				}
			}
		}

	}

	public static void main(String[] args) {
		char[][] matrix = new char[][] { { 'i', 'l', 'o', 'v', 'e' },
				{ 'd', 'i', 'n', 't', 'e' }, { 'n', 'e', 'w', 'e', 'p' },
				{ 'a', 'i', 'v', 'r', 'i' }, { 'm', 'a', 'x', 'e', 'c' } };
		MatrixSpiralRead matrixSpiralRead = new MatrixSpiralRead();
		matrixSpiralRead.spiralRead(matrix);

	}
}
