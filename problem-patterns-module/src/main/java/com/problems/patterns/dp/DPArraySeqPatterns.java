package com.problems.patterns.dp;

import java.util.Arrays;

/* All the below problems are solved in 4 approaches such as,
 *    1.Recursion -  Time: exponential time O(2^n), space complexity is O(n) which is used to store the recursion stack.
 *    2.DP: Top Down Approach or Memoization - Time and space complexity is same as memoization array size. Mostly same as Bottomup space and time.
 *    3.DP: Bottom Up Approach or Tabulation - Time & space complexity is O(n^2) or O(n)
 *    4.Memory Optimization - This will be same time as Bottom up approach, but space efficient. Eg: Two variable approach
 *    
 *  Imp Note: All the recursion problems needs below two points
 *  		1. Base case -> Its common for Recursion solution, memoization and bottom up approach
 *  		2. Recursive function(Sub problems solution)
 */
public class DPArraySeqPatterns {

	/***************************** Pattern 7: Array-Subsequence Probs *******************************/
	/* Longest Increasing Sequence:
	 * Eg: 
	 * 	Input: nums = [10,9,2,5,3,7,101,18], Output: 4
	 * 	Input: nums = [0,1,0,3,2,3], Output: 4
	 */
	/* Recursive Soln:
	 * 1. If the current element is greater than the previous element, then we can either pick it or dont pick it because we may get 
	 *    a smaller element somewhere ahead which is greater than previous and picking that would be optimal. So we try both options.
	 * 2. If the current element is smaller or equal to previous element, it can't be picked. 
	*/
	public int LIS1(int[] nums) {
		if (nums.length <= 1) return nums.length;
		//return lengthOfLIS(nums, 0, Integer.MIN_VALUE);
		return lengthOfLIS1(nums, 0, Integer.MIN_VALUE);
	}

	//Recursion1:
	public int lengthOfLIS1(int[] nums, int i, int prevNum) {
		if (i >= nums.length) return 0;

		// pick it if it is greater than previous picked element
		return Math.max(prevNum < nums[i] ? 1 + lengthOfLIS1(nums, i + 1, nums[i]) : 0,
				lengthOfLIS1(nums, i + 1, prevNum)); // return whichever choice gives max LIS
	}

	//Recursion2:
	public int lengthOfLIS(int[] nums, int i, int prevNum) {
		if (i >= nums.length) return 0;
		int taken = 0, notTaken = 0;
		if (prevNum < nums[i]) { // pick it if it is greater than previous picked element
			taken = 1 + lengthOfLIS(nums, i + 1, nums[i]);
		}
		notTaken = lengthOfLIS(nums, i + 1, prevNum);
		return Math.max(taken, notTaken); // return whichever choice gives max LIS
	}

	//Approach2: Top down approach 
	/*
	 * It wouldn't be scalable to store the state as (i, prev) because prev element can be any number in [-104, 104] 
	 * meaning we would need to declare a matrix dp[n][1e8] which won't be possibl
	 * DP with (i, prev) as state:
	 * Instead, we could store the state of (i, prev_i), where prev_i denotes the index of previous chosen element. 
	 * Thus we would use a dp matrix where dp[i][j] will denote the longest increasing subsequence from index i when
	 * previous chosen element's index is j.
	 */
	//TODO: This solution is not working. Revisit this
	public int LIS2(int[] nums) {
		if (nums.length <= 1) return nums.length;
		int[] dp = new int[nums.length + 1];
		Arrays.fill(dp, -1);
		return lengthOfLIS(nums, 0, -1, dp);
	}

	public int lengthOfLIS(int[] nums, int i, int prevIndex, int[] dp) {
		if (i >= nums.length) return 0;

		if (dp[prevIndex + 1] != -1) return dp[prevIndex + 1];

		int taken = 0, notTaken = 0;
		if (prevIndex == -1 || nums[prevIndex] < nums[i]) { // pick it if it is greater than previous picked element
			taken = 1 + lengthOfLIS(nums, i + 1, i);
		}
		notTaken = lengthOfLIS(nums, i + 1, prevIndex);
		return dp[prevIndex + 1] = Math.max(taken, notTaken); // return whichever choice gives max LIS
	}

	// Approach2: DP Approach : O(n^2)
	public int LIS3(int[] arr) {
		int n = arr.length;
		if (n <= 1) return n;
		int[] dp = new int[n];
		Arrays.fill(dp, 1);
		int max = dp[0];

		for (int i = 1; i < n; i++) {
			for (int j = 0; j < i; j++) {
				if (arr[j] < arr[i] && dp[i] < dp[j] + 1) {
					dp[i] = dp[j] + 1;
				}
			}
			max = Math.max(max, dp[i]);
		}
		return max;
	}

	// Binary Search Approach : O(nlogn)
	public int LIS4(int[] nums) {
		int[] dp = new int[nums.length];
		int size = 0;
		for (int num : nums) {
			//Here num is target, should find insert position for each num
			int index = searchInsert(dp, 0, size, num);

			// After finding the insert position store it in dp array
			dp[index] = num;

			//Size will be increased only insert position reaches the last index
			if (index == size) ++size;
		}
		return size;
	}

	//Here we should use BS insert position second approach, because here index h is arr.length   
	private int searchInsert(int[] nums, int l, int h, int target) {
		while (l < h) {
			int m = l + (h - l) / 2;
			if (target > nums[m]) {
				l = m + 1;
			} else {
				h = m;
			}
		}
		return l;
	}

	/* Longest Bitonic Subsequence:
	 * Given an array arr[0 … n-1] containing n positive integers, a subsequence of arr[] is called Bitonic if it is first increasing, then decreasing.
	 * Write a function that takes an array as argument and returns the length of the longest bitonic subsequence.
	 * Input arr[] = {1, 11, 2, 10, 4, 5, 2, 1}; Output: 6 (A Longest Bitonic Subsequence of length 6 is 1, 2, 10, 4, 2, 1)
	 * Input arr[] = {80, 60, 30, 40, 20, 10}; Output: 5 (A Longest Bitonic Subsequence of length 5 is 80, 60, 30, 20, 10)
	 */
	// Approach3: DP-Bottom Up Approach
	public int lbs3(int[] arr) {
		int n = arr.length;
		int[] lis = new int[n]; // Largest Increasing Sequence array
		int[] lds = new int[n];// Largest Decreasing Sequence array
		for (int i = 0; i < n; i++) {
			lis[i] = 1;
			lds[i] = 1;
		}

		// Largest Increasing Sequence logic
		for (int i = 1; i < n; i++)
			for (int j = 0; j < i; j++)
				if (arr[j] < arr[i] && lis[i] < lis[j] + 1) lis[i] = lis[j] + 1;

		// Largest Decreasing Sequence logic
		for (int i = n - 2; i >= 0; i--)
			for (int j = n - 1; j > i; j--)
				if (arr[j] < arr[i] && lds[i] < lds[j] + 1) lds[i] = lds[j] + 1;

		// Find the Bitonic value from LIS & LDS ( LIS+LDS-1)
		int max = Integer.MIN_VALUE;
		for (int i = 0; i < n; i++) {
			max = Math.max(max, lis[i] + lds[i] - 1);
		}
		return max;
	}

	/* Maximum Sum Increasing Subsequence:
	 * Given an array of n positive integers. Write a program to find the sum of maximum sum subsequence of the given array such that the
	 *  integers in the subsequence are sorted in increasing order. 
	 *  For example, if input is {1, 101, 2, 3, 100, 4, 5}, then output should be 106 (1 + 2 + 3 + 100),
	 *  			 if the input array is {3, 4, 5, 10}, then output should be 22 (3 + 4 + 5 + 10) 
	 *   			 if the input array is {10, 5, 4, 3}, then output should be 10
	 */
	// Approach3: DP-Bottom Up Approach
	public int MSIS3(int[] a) {
		int n = a.length;
		if (n <= 1) return n;
		int[] dp = new int[n], indexSeq = new int[n];
		int max = 0;
		for (int i = 0; i < n; i++) {
			dp[i] = a[i];
			indexSeq[i] = i;
		}
		for (int i = 1; i < n; i++) {
			for (int j = 0; j < i; j++) {
				if (a[j] < a[i] && dp[i] < dp[j] + a[i]) {
					dp[i] = dp[j] + a[i];
					indexSeq[i] = j;
				}
			}
			max = Math.max(max, dp[i]);
		}
		printMSIS(a, dp, indexSeq);
		return max;
	}

	private void printMSIS(int[] a, int[] msis, int[] indexSeq) {
		int max = Integer.MIN_VALUE, index = 0;
		for (int i = 0; i < msis.length; i++) {
			if (msis[i] > max) {
				max = msis[i];
				index = i; //Max Index
			}
		}
		int temp = max;
		while (temp > 0) {
			System.out.print(a[index] + " ");
			temp = temp - a[index];
			index = indexSeq[index];
		}
	}

}