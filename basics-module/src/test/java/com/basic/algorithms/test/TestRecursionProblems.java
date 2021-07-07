package com.basic.algorithms.test;

import java.util.Arrays;

import com.basic.algorithms.Recursion;

public class TestRecursionProblems extends Recursion {
	public static void main(String[] args) {
		TestRecursionProblems ob = new TestRecursionProblems();

		// Methods to understand the recursion concepts
		ob.testSingleRecursion();
		ob.testMultipleRecursion();
		// ob.testIterativeVsRecursive();
	}

	public void testSingleRecursion() {

		System.out.println("\nDecimal to Bin:" + decimalToBin(5));

		System.out.println("\nHead Recursive call:");
		increasingNumber(5);

		System.out.println("\nTail Recursive call:");
		decreasingNumber(5);
		// System.out.println("\nTotal no of calls in tailRecursiveCall: " + count);
	}

	public void testMultipleRecursion() {
		System.out.println("\nFib of value n:" + fibRecursive(6));
		// System.out.println("\nFib of value n:" + fibIterative(6));

		System.out.println("\nPrint Zig Zag:");
		printZigZag(3);

		int[] arr = { 0, 1, 89, 144, 2, 233, 5, 8, 377, 3, 1, 13, 21, 34, 55, 610 };
		System.out.println("Before sorting:");
		Arrays.stream(arr).forEach(a -> System.out.print(a + " "));
		System.out.println("\nAfter sorting:");
		Arrays.stream(arr).forEach(a -> System.out.print(a + " "));
	}

	public void testIterativeVsRecursive() {
		System.out.println("Single Iteration:");
		singleIteration(4);

		System.out.println("\nSingle Recursion:");
		singleRecursion(1, 4);

		System.out.println("\nDouble Iteration:");
		doubleIteration(5);

		System.out.println("\nDouble Recursion: ");
		doubleRecursion(4);
		System.out.println("No of method calls: " + dualRecursionCount);

		System.out.println("\nOne Iteration & Recursion:");
		oneIterAndRecursion(4);
		System.out.println("No of method calls: " + oneIterAndRecurCount);

		System.out.println("\nSum of all the double recursion combination: ");
		doubleRecursionWithSum(4);

		System.out.println("\nTriple Iteration:");
		tripleIteration(4);

		System.out.println("\nTriple Recursion: ");
		tripleRecursion(1, 4);
		System.out.println("\nNo of executions: " + tripleRecursionCount);
	}

}