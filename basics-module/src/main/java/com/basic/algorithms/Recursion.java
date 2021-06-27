package com.basic.algorithms;

import java.util.ArrayList;
import java.util.List;

/* 
 * Recursion:
 *  The process in which a function calls itself directly or indirectly is called recursion and the corresponding function is called as recursive
 *  function. Using recursive algorithm, certain problems can be solved quite easily. 
 *  Examples of such problems are Towers of Hanoi (TOH), Inorder/Preorder/Postorder Tree Traversals, DFS of Graph, etc.
 *   
 * Recursion Template
 *   - Base Case or Termination Condition
 *   - Any Logic: The problem can broken down into smaller problems of same type.
 *   - Recursive Function calls: It can be one or more time
 *   
 * Types:
 *  1.Linear or Single Recursion: A linear recursive function is a function that only makes a single call to itself each time the function runs
 *  	i. Head Recursion: First make a recursive call, then do the logic(any logic like add, div, sub, mul, print etc)
 *  	ii.Tail Recursion: First do a logic, then make recursive call.
 *  	Note: Head and Tail recursion applicable for all types of recursion such as linear, tree, mutual and nested.
 *  2.Tree or Multiple Recursion: If a recursive function calling itself for more than one time then it’s known as Tree Recursion.
 *  	- Functions with two recursive calls are referred to as binary recursive functions. No of function calls: 2^n
 *        Binary Recursion examples: Binary Tree traversals
 *      - Function calls itself more than twice. No of function calls: a^n
 *  4.Mutual or Indirect Recursion: In this recursion, there may be more than one functions and they are calling one another in a circular manner.
 *  5.Nested Recursion: In this recursion, a recursive function will pass the parameter as a recursive call.That means “recursion inside recursion”. 
 *   
 * When you hear a problem beginning with the following statements, it's often (though not always) a good candidate for recursion:
 *  - "Design an algorithm to compute the nth .. :;  
 *  - "Write code to list the first n .. :; 
 *  - "Implement a method to compute all..:; and so on.
 *  
 * Recursive Problems Categorized by "ByteByBye" website:
 *  - Iteration - Iterate over a variety of data structures using recursion, both in one and multiple dimensions; 
 *  	Eg: Insert Element at the Bottom of a Stack, Generating All Substrings of a String, Flattening a 2D Array, and more
 *  - SubProblems - Most fundamental pattern in all of recursion: Subproblems; 
 *  	Eg: Stair Stepping, Towers of Hanoi, Is String a Palindrome, and more
 *  - Selection - Pattern related to DP & Combination Problems. 
 *  	Eg: Find All Combinations, 0-1 Knapsack, String Interleaving, and more
 *  - Ordering - Pattern related to permutations & similar problems;  
 *  	Eg: Find All Permutations, N-digit Numbers, BST Arrays, and more
 *  - Divide & Conquer - solve problems by breaking them into smaller pieces, grouping subpattern of difficult problems; 
 *  	Eg: Binary Search, Unique BSTs, String Compression, Rotated Arrays, and more
 *  - Depth First Search - Mostly applied to Tree, Graph & Matrix; 
 *  	Eg: DFS in Trees and Graphs, Find all Combinations via DFS, and more
 
 */
public class Recursion {

	MathProblems math = new MathProblems();

	SortingAlgorithms sort = new SortingAlgorithms();

	public static int count = 0;

	/**************** Methods to understand the recursion concepts *******************/
	public void recursionMethodConsolidations(int n, int[] a) {
		//1.Linear Recursion
		//Head Vs Tail recursion types
		increasingNumber(n);
		decreasingNumber(n);
		//More examples in Linear Recursion
		factorial(n);
		sqrt(9, 2);
		gcd(10, 5);
		decimalToBin(n);

		//2.Tree or Multiple Recursion
		fibRecursive(n);
		quickSort(a);
		printZigZag(n);
	}

	/*************** 1.Single or Linear Recursion ********************/
	//Head recursive type
	public void increasingNumber(int n) {
		//Base case
		if (n == 0) return;

		//Recursive Function call
		increasingNumber(n - 1);

		//Any logic: here just printing the number 
		System.out.print(n + " "); //Result: 1 2 3 4 5
	}

	// Tail recursive type
	public void decreasingNumber(int n) {
		//Base case
		if (n == 0) return;

		//Any logic: here just printing the number 
		System.out.print(n + " "); //Result: 5 4 3 2 1

		//Recursive Function call
		decreasingNumber(n - 1);
	}

	//Example for linear recursion: Factorial , sqrt, gcd, decimalToBinary etc
	public void factorial(int n) {
		math.factorial1(n);
	}

	//Compute the square root of a number using Newton's method (assume EPSILON to be a very small number close to 0):
	public double sqrt(double x, double a) {
		return math.mySqrt4(x, a);
	}

	//Example for a tail recursive function is GCD
	public void gcd(int a, int b) {
		math.gcd1(a, b);
	}

	// Example2: Decimal to Binary Conversion
	public int decimalToBin(int n) {
		return math.decimalToBin2(n);
	}

	/*************** 2.Multiple or Tree Recursion ********************/
	// Example1: Fibonacci Series; Time-O(2^n)
	public long fibRecursive(int n) {
		if (n <= 1) return n;
		return fibRecursive(n - 1) + fibRecursive(n - 2);
	}

	// Fibonacci series: 0,1,1,2,3,5,8,13,21,34...
	//Time: O(n)
	public long fibIterative(int n) {
		if (n <= 1) return n;
		long f1 = 0, f2 = 1;
		System.out.print("Fibonacci series: ");
		System.out.print(f1 + " " + f2 + " ");
		int i = 2;
		long result = 0;
		while (i++ <= n) {
			result = f1 + f2;
			f1 = f2;
			f2 = result;
			System.out.print(result + " ");
		}
		return result;
	}

	// Example2: Quick Sort
	// Binary Recursion: Tail recursion type
	public void quickSort(int[] a) {
		sort.quickSort(a);
	}

	// Print Zig Zag
	public void printZigZag(int n) {
		if (n == 0) return;

		System.out.println("Pre: " + n);

		printZigZag(n - 1);

		System.out.println("In: " + n);

		printZigZag(n - 1);

		System.out.println("Post: " + n);
	}

	//TODO:  Revisit and delete below

	/*************** 4.Iterative Vs Recursive ********************/
	// Single Iteration
	public void singleIteration(int n) {
		for (int i = 1; i <= n; i++)
			System.out.print(i + ", ");
	}

	// Single Recursion
	public void singleRecursion(int i, int n) {
		if (i > n) return;
		System.out.print(i + ", ");
		singleRecursion(i + 1, n);
	}

	/*Double Iteration -> No of executions: n(n+1)/2 for this case. It varies based on the condition. 
	 * By default, no of execution: n*n times*/
	public void doubleIteration(int n) {
		for (int i = 1; i <= n; i++)
			for (int j = i; j <= n; j++)
				System.out.print(j + ", ");
	}

	public static int dualRecursionCount = 0;

	/*Two Recursion -> No of executions: (2^n) - 1 or PrevValue+2^(n-1); Series: 1,3,7,15,31,63....*/
	public void doubleRecursion(int n) {
		doubleRecursion(0, n);
	}

	private void doubleRecursion(int index, int n) {
		if (index >= n) return;

		dualRecursionCount++;
		System.out.print(index + ", ");
		doubleRecursion(index + 1, n);
		// System.out.println();
		doubleRecursion(index + 1, n);
	}

	public static int oneIterAndRecurCount = 0;

	// One iteration & Recursion: Same result like double recursion method
	public void oneIterAndRecursion(int n) {
		// oneIterAndRecursion(0, n);
		oneIterAndRecursion2(n);
	}

	public void oneIterAndRecursion(int index, int n) {
		if (index >= n) return;
		// System.out.print(i + ", ");
		for (int j = index; j < n; j++) {
			oneIterAndRecurCount++;
			System.out.print(j + ", ");
			oneIterAndRecursion(j + 1, n);
		}
	}

	public void oneIterAndRecursion2(int n) {
		if (n <= 0) return;
		// System.out.print(i + ", ");
		for (int i = 0; i < n; i++) {
			oneIterAndRecurCount++;
			System.out.print(i + ", ");
			oneIterAndRecursion2(n - i - 1);
		}
	}

	// Example for double recursion: Sum of n combinations
	public void doubleRecursionWithSum(int n) {
		System.out.println("Sum: ");
		// sumOfSequences(1, 4, 0);
		List<Integer> seq = new ArrayList<>();
		sumOfSequences(1, 4, 0, seq);
	}

	public void sumOfSequences(int i, int n, int sum) {
		if (i > n) {
			System.out.print(sum + ", ");
			return;
		}

		sum += i;
		// System.out.print(i + "-" + sum + ", ");
		sumOfSequences(i + 1, n, sum);
		sumOfSequences(i + 1, n, sum - i);
	}

	public void sumOfSequences(int i, int n, int sum, List<Integer> seq) {
		if (i > n) {
			seq.stream().forEach(k -> System.out.print(k + ","));
			System.out.print("=" + sum + "; ");
			return;
		}

		sum += i;
		seq.add(i);
		// System.out.print(i + "-" + sum + ", ");
		sumOfSequences(i + 1, n, sum, seq);

		if (!seq.isEmpty()) seq.remove(seq.size() - 1);

		sumOfSequences(i + 1, n, sum - i, seq);
	}

	// Triple Iteration
	public void tripleIteration(int n) {
		for (int i = 1; i <= n; i++)
			for (int j = i; j <= n; j++)
				for (int k = j; k <= n; k++)
					System.out.print(k + ", ");
	}

	public static int tripleRecursionCount = 0;

	/*Three Recursion -> No of executions: PrevValue+3^(n-1); Series: 1,4,13,40,121,364....*/
	public void tripleRecursion(int i, int n) {
		if (i > n) return;

		System.out.print(i + ", ");
		tripleRecursionCount++;
		tripleRecursion(i + 1, n);
		tripleRecursion(i + 1, n);
		tripleRecursion(i + 1, n);
	}
}