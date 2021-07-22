package com.basic.datastructures;

//TODO: Think about what to be here????
public class Matrix {

	//2D Array basics
	public void basicApis() {
		//2D array initialization
		//2D array sorting 
		//2D array swapping 
	}

	/*************** Matrix Basic Problems ************************/
	public int[][] transpose(int[][] A) {
		int row = A.length, col = A[0].length;
		int[][] transpose = new int[col][row];
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				transpose[j][i] = A[i][j];
			}
		}
		return transpose;
	}

	// Matrix Multiplication
	public int[][] matrixMul(int[][] a, int[][] b) {
		int m1 = a.length, n1 = a[0].length;
		int m2 = b.length, n2 = b[0].length;
		int[][] result = new int[m1][n2];

		/* Matrix Mul condition: 
		   		Given Matrix: a(m1 X n1), b(m2 X n2)
		   		To perform multiplication, n1 == m2  and 
		   		Result size: m1 X n2
		 */
		if (n1 != m2) return result;

		for (int i = 0; i < m1; i++)
			for (int j = 0; j < n2; j++)
				for (int k = 0; k < n1; k++) // k -> n1 or m2
					result[i][j] += a[i][k] * b[k][j];

		return result;
	}

	/*Sparse Matrix Multiplication:
	 * Optimized Method: We can see that when a_ik is 0, there is no need to compute b_kj. So we switch the inner two 
	 * loops and add a 0-checking condition. 
	 * Since the matrix is sparse, time complexity is ~O(n^2) which is much faster than O(n^3).
	 */
	public int[][] sparseMatrixMul(int[][] a, int[][] b) {
		/*int m = a.length, n = b[0].length;
		int[][] result = new int[m][n];*/

		int m1 = a.length, n1 = a[0].length;
		int m2 = b.length, n2 = b[0].length;
		int[][] result = new int[m1][n2];

		if (n1 != m2) return result;

		for (int i = 0; i < m1; i++) {
			for (int k = 0; k < n1; k++) { // k -> n1 or m2
				if (a[i][k] == 0) continue;
				for (int j = 0; j < n2; j++)
					result[i][j] += a[i][k] * b[k][j];
			}
		}

		return result;
	}

	/*************** Matrix Traversals or Print ************************/

	/* Spiral Matrix:
	 * Given a matrix of m x n elements (m rows, n columns), return all elements of the matrix in spiral order.
	 * 	Example 1:
		Input:
			[
			[ 1, 2, 3 ],
			[ 4, 5, 6 ],
			[ 7, 8, 9 ]
			]
		Output: [1,2,3,6,9,8,7,4,5]
	 */

	public int[] spiralOrder(int[][] matrix) {
		int r = matrix.length, c = matrix[0].length;
		int[] result = new int[r * c];

		if (matrix.length == 0 || matrix[0].length == 0) return result;
		int left = 0, right = c - 1, top = 0, bottom = r - 1, index = 0;

		while (top <= bottom && left <= right) {
			for (int j = left; j <= right; j++)
				result[index++] = matrix[top][j];
			top++;

			for (int i = top; i <= bottom; i++)
				result[index++] = matrix[i][right];
			right--;

			if (top > bottom || left > right) break;

			for (int j = right; j >= left; j--)
				result[index++] = matrix[bottom][j];
			bottom--;

			for (int i = bottom; i >= top; i--)
				result[index++] = matrix[i][left];
			left++;
		}

		return result;
	}

	/*
	 * Spiral Matrix II: 
	 * Given a positive integer n, generate a square matrix filled with elements from 1 to n2 in spiral order.
	 */
	public int[][] generateMatrix(int n) {
		int l = 0, r = n - 1, t = 0, b = n - 1;
		int[][] mat = new int[n][n];
		int val = 1;
		while (l <= r && t <= b) {
			for (int j = l; j <= r; j++)
				mat[t][j] = val++;
			t++;
			for (int i = t; i <= b; i++)
				mat[i][r] = val++;
			r--;
			if (l > r || t > b) break;
			for (int j = r; j >= l; j--)
				mat[b][j] = val++;
			b--;
			for (int i = b; i >= t; i--)
				mat[i][l] = val++;
			l++;
		}

		return mat;
	}

	/* Diagonal Traverse/Print a 2D array in Diagonal ZigZag order
	 * Solution:
	 * 	- Regular Case: Traversal should be up or down diagonally based on sum of indices(even or odd)
	 * 		- Diagonally upward: r--, c++
	 *      - Diagonally downward: r++, c--;
	 *  - Corner Case: Whenever we reach an array boundary(r=0 or m-1 & c=0 or m-1), we need to do two things:
	 *  	1.whenever we reach a row boundary (first or last row) , we shift our path right.
	 *   		- if r==0 or r==m-1, then move right, i.e c++
	 *   	2.whenever we reach a column boundary (first or last column) , we shift our path down.
	 *   		- if c==0 or c==n-1, then move down, i.e r--
	 */
	public int[] findDiagonalOrder(int[][] matrix) {
		if (matrix == null || matrix.length == 0) return new int[0];

		int r = 0, c = 0;
		int m = matrix.length, n = matrix[0].length;

		int[] result = new int[m * n];

		for (int i = 0; i < result.length; i++) {
			result[i] = matrix[r][c];
			//If sum of indices are even, then should traverse upward diagonally.
			if ((r + c) % 2 == 0) {
				if (c == n - 1) { //Move down
					r++;
				} else if (r == 0) { // Move right
					c++;
				} else { //upward diagonal traversal
					r--;
					c++;
				}
			} else {//If sum of indices are odd, then should traverse downward diagonally.
				if (r == m - 1) { //Move right
					c++;
				} else if (c == 0) { //Move down
					r++;
				} else { //downward diagonal traversal
					r++;
					c--;
				}
			}
		}

		return result;
	}

	// Diagonal Traverse - Only upward
	//Eg: https://www.geeksforgeeks.org/zigzag-or-diagonal-traversal-of-matrix/
	public void diagonalUpwarddirection(int[][] matrix) {
		int r = matrix.length, c = matrix[0].length;

		for (int k = 0; k < r; k++) {
			int i = k, j = 0;
			while (i >= 0 && j <= k) {
				System.out.print(matrix[i][j] + " ");
				i--;
				j++;
			}
			System.out.println();
		}

		for (int k = 1; k < c; k++) {
			int i = r - 1, j = k;
			while (i >= 0 && j < c) {
				System.out.print(matrix[i][j] + " ");
				i--;
				j++;
			}
			System.out.println();
		}
	}

}