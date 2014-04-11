package Recursion;

import java.util.ArrayList;

/* 
 * 8.2
 * Imagine a robot sitting on the upper left hand corner of an NxN grid The robot can only move in two directions: right and down 
 * How many possible paths are there for the robot?
 */

public class RobotPaths {
	int size;
	int[][] map;
	ArrayList<int[]> trace = new ArrayList<int[]>();
	ArrayList<ArrayList<int[]>> allTraces = new ArrayList<ArrayList<int[]>>();

	public RobotPaths(int n) {
		this.size = n;
		this.map = new int[n][n];
		for (int row = 0; row < size; row++) {
			for (int column = 0; column < size; column++) {
				map[row][column] = 1;
			}
		}
		// set hindrances
		// map[row][column] = 0;
	}

	public boolean searchPaths(int x, int y) {
		int[] step = new int[] { x, y };
		trace.add(step);
		if (x == 0 && y == 0) {
			ArrayList<int[]> traceCopy = new ArrayList<int[]>();
			traceCopy.addAll(trace);
			allTraces.add(traceCopy);
			// trace.clear();
			return true;
		}
		boolean success = false;
		// back trace right
		if (x > 0 && isClear(x - 1, y)) {
			success = searchPaths(x - 1, y);
			int[] backStep = new int[] { x - 1, y };
			// get rid of steps in previous path before searching a new path
			trace = removeStep(trace, backStep);
		}
		// back trace down
		if (/* !success && */y > 0 && isClear(x, y - 1)) {
			success = searchPaths(x, y - 1);
			int[] backStep = new int[] { x, y - 1 };
			trace = removeStep(trace, backStep);
		}
		// dead end, wrong step
		if (!success) {
			trace.remove(step);
		}
		return success;
	}

	public boolean isClear(int x, int y) {
		if (map[x][y] == 1) {
			return true;
		}
		return false;
	}

	public ArrayList<int[]> removeStep(ArrayList<int[]> trace, int[] discardStep) {
		ArrayList<int[]> newTrace = new ArrayList<int[]>();
		for (int[] oneStep : trace) {
			if (!java.util.Arrays.equals(oneStep, discardStep)) {
				newTrace.add(oneStep);
			}
		}
		return newTrace;
	}

	public void showTrace() {
		System.out.println("# of paths: " + allTraces.size());
		System.out.println("-------------------------------");
		for (ArrayList<int[]> oneTrace : allTraces) {
			for (int i = oneTrace.size() - 1; i > 0; i--) {
				int j = i - 1;
				if (oneTrace.get(j)[0] - oneTrace.get(i)[0] == 1) {
					System.out.print("(" + oneTrace.get(i)[0] + ","
							+ oneTrace.get(i)[1] + ")" + "→ ");
				} else if (oneTrace.get(j)[1] - oneTrace.get(i)[1] == 1) {
					System.out.print("(" + oneTrace.get(i)[0] + ","
							+ oneTrace.get(i)[1] + ")" + "↓ ");
				}
			}
			System.out.print("(" + oneTrace.get(0)[0] + ","
					+ oneTrace.get(0)[1] + ")");
			System.out.println();
		}
	}

	public static void main(String[] args) {
		int n = 5;
		RobotPaths robotPaths = new RobotPaths(n);
		robotPaths.searchPaths(n - 1, n - 1);
		robotPaths.showTrace();

	}

}
