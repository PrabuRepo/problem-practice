package com.consolidated.problems.algorithms;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/*
 * Backtracking Problems
 */
public class Backtracking {
	/********************* Backtracking Patterns – Auxiliary Buffer ***********************/

	/********************* Backtracking Other Problems ***********************/
	// Flip Game I, II:
	/*
	 * Flip Game:  You are playing the following Flip Game with your friend: Given a string that contains only these two
	 * characters: + and -, you and your friend take turns to flip two consecutive "++" into "--". The game ends when a
	 * person can no longer make a move and therefore the other person will be the winner. Write a function to compute
	 * all possible states of the string after one valid move. 
	 * For example, given s = "++++", after one move, it may become one of the following states: 
	 * [ "--++", "+--+", "++--" ] 
	 * If there is no valid move, return an empty list [].
	 */
	// Flip Game I
	public List<String> generatePossibleNextMoves1(String s) {
		List<String> res = new ArrayList<>();
		int n = s.length();
		for (int i = 1; i < n; i++) {
			if (s.charAt(i - 1) == '+' && s.charAt(i) == '+') {
				res.add(s.substring(0, i - 1) + "--" + s.substring(i + 1, n));
			}
		}
		return res;
	}

	public List<String> generatePossibleNextMoves2(String s) {
		List<String> res = new ArrayList<String>();
		if (s == null) return res;
		char[] arr = s.toCharArray();
		for (int i = 0; i < arr.length - 1; i++) {
			if (arr[i] == arr[i + 1] && arr[i] == '+') {
				arr[i] = '-';
				arr[i + 1] = '-';
				res.add(new String(arr));
				arr[i] = '+';
				arr[i + 1] = '+';
			}
		}
		return res;
	}

	// Flip Game II:
	// Approach1: Backtracking Solution
	public boolean canWin1(String s) {
		if (s == null || s.length() < 2) return false;
		for (int i = 0; i < s.length() - 1; i++) {
			if (s.startsWith("++", i)) {
				String subStr = s.substring(0, i) + "--" + s.substring(i + 2);
				if (!canWin1(subStr)) return true;
			}
		}
		return false;
	}

	// Approach2: Optimization: DP+ memory search
	public boolean canWin2(String s) {
		if (s == null || s.length() < 2) return false;
		Map<String, Boolean> map = new HashMap<>();
		return helper(s, map);
	}

	public boolean helper(String s, Map<String, Boolean> map) {
		if (map.containsKey(s)) return map.get(s);
		for (int i = 0; i < s.length() - 1; i++) {
			if (s.startsWith("++", i)) {
				String subStr = s.substring(0, i) + "--" + s.substring(i + 2);
				if (!helper(subStr, map)) {
					map.put(s, true);
					return true;
				}
			}
		}
		map.put(s, false);
		return false;
	}

	// Minimum Unique Word Abbreviation - Heap/Trie/Bactracking

	// N-Queens, N-Queens II/Eight Queens
	// N Queen Problem
	public boolean solveNQ() {
		int board[][] = { { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 } };

		if (solveNQUtil(board, 0) == false) {
			System.out.print("Solution does not exist");
			return false;
		}

		printSolution(board);
		return true;
	}

	// A recursive utility function to solve N Queen problem
	boolean solveNQUtil(int board[][], int col) {
		int N = board.length;
		// base case: If all queens are placed then return true
		if (col >= N) return true;

		// Consider this column and try placing this queen in all rows one by one
		for (int i = 0; i < N; i++) {
			// Check if the queen can be placed on board[i][col]
			if (isSafeBoard(board, i, col)) {
				board[i][col] = 1; // Place this queen in board[i][col]

				if (solveNQUtil(board, col + 1) == true) // recur to place rest of the queens
					return true;

				// If placing queen in board[i][col] doesn't lead to a solution then remove queen from board[i][col]
				board[i][col] = 0; // BACKTRACK
			}
		}
		return false;
	}

	/* A utility function to check if a queen can be placed on board[row][col]. Note that this function is called when "col" queens
	 * are already placed in columns from 0 to col -1. So we need to check only left side for attacking queens */
	boolean isSafeBoard(int board[][], int row, int col) {
		int N = board.length;
		int i, j;

		/* Check this row on left side */
		for (i = 0; i < col; i++)
			if (board[row][i] == 1) return false;

		/* Check upper diagonal on left side */
		for (i = row, j = col; i >= 0 && j >= 0; i--, j--)
			if (board[i][j] == 1) return false;

		/* Check lower diagonal on left side */
		for (i = row, j = col; j >= 0 && i < N; i++, j--)
			if (board[i][j] == 1) return false;

		return true;
	}

	void printNQueenSolution(int board[][]) {
		int N = board.length;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++)
				System.out.print(" " + board[i][j] + " ");
			System.out.println();
		}
	}

	// Solve the Sudoku/Sudoku Solver
	public boolean solveSudoku(int[][] grid) {
		int n = grid.length; // Row & Col size is same

		// Find the unassigned Location
		int row = -1, col = -1;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				if (grid[i][j] == 0) {
					row = i;
					col = j;
				}

		// All the boxes are filled; success
		if (row == -1 && col == -1) return true;

		for (int i = 1; i <= n; i++) {
			if (isSafe(grid, row, col, i)) {
				grid[row][col] = i;

				if (solveSudoku(grid)) // Success
					return true;

				// Backtracking: Failure, reset the value and try again
				grid[row][col] = 0;
			}
		}

		return false;
	}

	public boolean isAnyUnassignedLocation(int[][] grid) {
		int n = grid.length;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				if (grid[i][j] == 0) return true;
		return false;
	}

	public boolean isSafe(int[][] grid, int row, int col, int num) {
		int boxSize = 3;
		if (!usedInRow(grid, row, num) && !usedInCol(grid, col, num)
				&& !usedInBox(grid, row - row % boxSize, col - col % boxSize, num))
			return true;
		return false;
	}

	public boolean usedInRow(int[][] grid, int row, int num) {
		for (int j = 0; j < grid.length; j++)
			if (grid[row][j] == num) return true;
		return false;
	}

	public boolean usedInCol(int[][] grid, int col, int num) {
		for (int i = 0; i < grid.length; i++)
			if (grid[i][col] == num) return true;
		return false;
	}

	public boolean usedInBox(int[][] grid, int boxStartRow, int boxStartCol, int num) {
		int boxSize = 3;
		for (int i = 0; i < boxSize; i++)
			for (int j = 0; j < boxSize; j++)
				if (grid[i + boxStartRow][j + boxStartCol] == num) return true;
		return false;
	}

	public void printGrid(int[][] grid) {
		int N = grid.length;
		for (int row = 0; row < N; row++) {
			for (int col = 0; col < N; col++)
				System.out.print(grid[row][col] + "  ");
			System.out.println();
		}
	}

	// Crossword Puzzle
	static final int SIZE = 10;
	static final int[] R_OFFSETS = { 0, 1 };
	static final int[] C_OFFSETS = { 1, 0 };

	// Crossword Puzzle: DFS & Backtracking
	public char[][] solvePuzzle(char[][] grid, String words) {
		return search(grid, Arrays.stream(words.split(";")).collect(Collectors.toSet()), 0, 0, 0);
	}

	public char[][] search(char[][] grid, Set<String> words, int r, int c, int direction) {
		if (r == SIZE) return grid;
		if (c == SIZE) return search(grid, words, r + 1, 0, 0);
		if (direction == R_OFFSETS.length) return search(grid, words, r, c + 1, 0);

		// Count the length of the path in the grid
		int insertLength = countInsertLength(grid, r, c, direction);

		if (insertLength > 1) {
			for (String word : new ArrayList<>(words)) {
				// Validate the word can be inserted in grid
				if (canInsertWord1(grid, r, c, direction, insertLength, word)) {
					List<Integer> insertOffsets = new ArrayList<Integer>();

					for (int i = 0; i < insertLength; i++) {
						int row = r + R_OFFSETS[direction] * i;
						int col = c + C_OFFSETS[direction] * i;

						if (grid[row][col] == '-') {
							grid[row][col] = word.charAt(i);

							insertOffsets.add(i);
						}
					}
					words.remove(word);

					char[][] subResult = search(grid, words, r, c, direction + 1);

					if (subResult != null) return subResult;

					// Backtracking: Reassign the values
					words.add(word);

					for (int insertOffset : insertOffsets) {
						// Calculate row & col using prev offset
						int row = r + R_OFFSETS[direction] * insertOffset;
						int col = c + C_OFFSETS[direction] * insertOffset;

						grid[row][col] = '-';
					}
				}
			}

			return null;
		} else {
			return search(grid, words, r, c, direction + 1);
		}
	}

	public int countInsertLength(char[][] grid, int r, int c, int direction) {
		int prevRow = r - R_OFFSETS[direction];
		int prevCol = c - C_OFFSETS[direction];

		if (prevRow >= 0 && prevRow < SIZE && prevCol >= 0 && prevCol < SIZE && grid[prevRow][prevCol] != '+') return 0;

		int insertLength = 0;
		while (r >= 0 && r < SIZE && c >= 0 && c < SIZE && grid[r][c] != '+') {
			insertLength++;
			r += R_OFFSETS[direction];
			c += C_OFFSETS[direction];
		}
		return insertLength;
	}

	public boolean canInsertWord1(char[][] grid, int r, int c, int direction, int insertLength, String word) {
		if (word.length() != insertLength) return false;

		for (int k = 0; k < word.length(); k++) {
			int row = r + R_OFFSETS[direction] * k;
			int col = c + C_OFFSETS[direction] * k;
			if (grid[row][col] != '-' && grid[row][col] != word.charAt(k)) return false;
		}

		return true;
	}

	public boolean canInsertWord2(char[][] grid, int r, int c, int direction, int insertLength, String word) {

		return word.length() == insertLength && IntStream.range(0, word.length()).allMatch(k -> {
			int row = r + R_OFFSETS[direction] * k;
			int col = c + C_OFFSETS[direction] * k;

			return grid[row][col] == '-' || grid[row][col] == word.charAt(k);
		});
	}

	// The Knight’s tour problem
	public void knightTour(int N) {
		int moveCount = 1;
		int[][] table = new int[N][N];

		int[] xMove = { 2, 1, -1, -2, -2, -1, 1, 2 };
		int[] yMove = { 1, 2, 2, 1, -1, -2, -2, -1 };

		/*Why this possibilities are not working?
		 * int[] xMove = { 2, 2, -2, -2, 1, 1, -1, -1 };
		int[] yMove = { 1, -1, 1, -1, 2, -2, 2, -2 };*/

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				table[i][j] = -1;

		table[0][0] = 0;

		knightTourUtil(0, 0, N, moveCount, xMove, yMove, table);

		printSolution(table, N);
	}

	static int count = 0;

	private boolean knightTourUtil(int x, int y, int N, int moveCount, int[] xMove, int[] yMove, int[][] table) {
		int nextX, nextY;
		count++;
		if (moveCount == (N * N)) return true;

		for (int i = 0; i < N; i++) {
			nextX = x + xMove[i];
			nextY = y + yMove[i];
			if (isSafe(nextX, nextY, N, table)) {
				// System.out.println("if-> X:" + nextX + " Y:" + nextY);
				table[nextX][nextY] = moveCount;
				if (knightTourUtil(nextX, nextY, N, moveCount + 1, xMove, yMove, table)) {
					return true;
				} else {
					// System.out.println("else-> X:" + nextX + " Y:" + nextY);
					table[nextX][nextY] = -1;
				}
			}
		}

		return false;
	}

	private boolean isSafe(int nextX, int nextY, int N, int[][] table) {
		return (nextX >= 0 && nextX < N && nextY >= 0 && nextY < N && table[nextX][nextY] == -1);
	}

	private void printSolution(int[][] table, int N) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				System.out.print(table[i][j] + "  ");
			}
			System.out.println();
		}
	}

	// Rat in a Maze
	public boolean solveMaze(int maze[][]) {
		int sol[][] = { { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 } };

		if (solveMazeUtil(maze, 0, 0, sol) == false) {
			System.out.print("Solution doesn't exist");
			return false;
		}

		printSolution(sol);
		return true;
	}

	// A recursive utility function to solve Maze problem
	boolean solveMazeUtil(int maze[][], int x, int y, int sol[][]) {
		int N = maze.length;
		if (x == N - 1 && y == N - 1) { // Base case: When it reaches the end of row & col, return true
			sol[x][y] = 1;
			return true;
		}

		// Check if maze[x][y] is valid
		if (isSafe(maze, x, y) == true) {
			sol[x][y] = 1; // mark x,y as part of solution path

			if (solveMazeUtil(maze, x + 1, y, sol)) // Move forward in x direction
				return true;

			if (solveMazeUtil(maze, x, y + 1, sol)) // Move down in y direction
				return true;

			sol[x][y] = 0; // If none of the above movements works then BACKTRACK: un mark x,y as part of solution path
			return false;
		}

		return false;
	}

	boolean isSafe(int maze[][], int x, int y) {
		int N = maze.length;
		// if (x,y outside maze) return false
		return (x >= 0 && x < N && y >= 0 && y < N && maze[x][y] == 1);
	}

	public void printSolution(int sol[][]) {
		int N = sol.length;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++)
				System.out.print(" " + sol[i][j] + " ");
			System.out.println();
		}
	}

}