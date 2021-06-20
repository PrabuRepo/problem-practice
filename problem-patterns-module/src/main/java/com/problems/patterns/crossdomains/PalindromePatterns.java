package com.problems.patterns.crossdomains;

import com.basic.algorithms.MathProblems;
import com.common.model.ListNode;
import com.problems.patterns.BacktrackingPatterns;
import com.problems.patterns.StringProblems;
import com.problems.patterns.dp.DPStringPatterns;
import com.problems.patterns.ds.FastAndSlowPtrPatterns;

/*
 * Palindrome Problem Patterns have three catagories,
 * 	1. Number: Here reverse the whole number and compare
 * 	2. LinkedList: Here reverse the second half and compare with first half of the list
 *  3. String: Compare char from both side and move towards mid to check the palindrome string
 */
public class PalindromePatterns {

	private StringProblems stringProblems;

	private MathProblems mathProblems;

	private FastAndSlowPtrPatterns fastAndSlowPtrPatterns;

	private DPStringPatterns dpStringPatterns;

	private BacktrackingPatterns backtrackingPatterns;

	/************************** Palindrome Number Problems ***********************************/
	//Palindrome Number  
	public void palindromeNumber(int num) {
		mathProblems.isPalindrome(num);
	}

	/* Find the Closest Palindrome
	 *  Example: 
	 *  	Input: 7599; Output: 7557
	 *  	Input: 7501; Output: 7447
	 *  	Input: 9999; Output: 10001
	 *  	Input: 10000; Output: 9999
	 */
	public String closestPalindrome(String n) {
		Long num = Long.valueOf(new String(n));

		//Order used to eliminate second half of digits
		int order = (int) Math.pow(10, n.length() / 2);

		/* 3 cases to mid: Example for 7599,
		 *   1. Same mid: 7557
		 *   2. Increase one to mid to get 7667
		 *   3. Decrease one to mid to get 7447
		 */
		//1.Same mid
		Long mirrorNum = mirror(num);

		//2.Increase 1 to mid and also handles input like 9,99,999...
		Long mirrorLarger = mirror((num / order) * order + order);

		//3.Decrease 1 to mid and also handles input like 10,100,1000...
		Long mirrorSmaller = mirror((num / order) * order - 1);

		//Below logic to find the closest to given number
		if (mirrorNum > num) {
			mirrorLarger = (long) Math.min(mirrorLarger, mirrorNum);
		} else if (mirrorNum < num) {
			mirrorSmaller = (long) Math.max(mirrorSmaller, mirrorNum);
		}

		Long closestValue = Math.abs(num - mirrorSmaller) <= Math.abs(mirrorLarger - num) ? mirrorSmaller
				: mirrorLarger;
		return String.valueOf(closestValue);
	}

	/*
	 * Example for Mirror:
	 * 	7599 become 7557
	 *  7400 become 7447
	 *  7600 become 7669
	 */
	private Long mirror(Long ans) {
		char[] a = String.valueOf(ans).toCharArray();
		int i = 0;
		int j = a.length - 1;
		while (i < j) {
			a[j--] = a[i++];
		}
		return Long.valueOf(new String(a));
	}

	/**************************** Palindrome Linked List Problems *********************************/
	//Palindrome Linked List
	public void palindromeLinkedList(ListNode head) {
		fastAndSlowPtrPatterns.isPalindrome1(head);
		fastAndSlowPtrPatterns.isPalindrome2(head);
	}

	/***************************** Palindrome String Problems ********************************/
	//Valid Palindrome I
	public void validPalindromeI(String str) {
		stringProblems.isPalindrome11(str);
		stringProblems.isPalindrome12(str);
		stringProblems.isPalindrome13(str);
	}

	//Valid Palindrome II
	public void validPalindromeII(String str) {
		stringProblems.validPalindrome2(str);
	}

	//Longest Palindromic Subsequence   
	public void palindromicSubsequence(String str) {
		dpStringPatterns.lps1(str);
		dpStringPatterns.lps3(str);
	}

	//Longest Palindromic Substring   
	public void palindromicSubstring(String str) {
		dpStringPatterns.lpSubstr1(str);
		dpStringPatterns.lpSubstr3(str);
	}

	//Palindrome Partitioning
	public void palindromePartitioning(String str) {
		backtrackingPatterns.partition1(str);
		backtrackingPatterns.partition2(str);
	}

	//Palindrome Partitioning II
	public void palindromePartitioningII(String str) {
		dpStringPatterns.minCut1(str);
		dpStringPatterns.minCut2(str);
	}

	//These problems are modification above LPS pattern
	//Minimum Deletions in a String to make it a Palindrome
	//Form a Palindrome(min no of chars needed to form palindrome)
	//Count of Palindromic Subsequence
	//Count of Palindromic Substrings/Palindromic Substrings

	//TODO: Solve below problems
	//Trie: Palindrome Pairs    
	//KMP Algorithm: Shortest Palindrome 

}