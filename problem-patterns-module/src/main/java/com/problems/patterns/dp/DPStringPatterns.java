package com.problems.patterns.dp;

import java.util.Arrays;

import com.common.utilities.Utils;
import com.problems.patterns.BacktrackingPatterns;

/*
 * Patterns covered in this class are,
 * 	1.Palindromic substring/subseq Probs
 *  2.Substring/Subsequence Probs - LCSs
 */
public class DPStringPatterns {

	private BacktrackingPatterns backtrackingPatterns;

	/********************* Pattern 5: String-Palindromic substring/subseq Probs ***********************/
	// Longest Palindromic Subsequence:
	// 1.Recursion Approach
	public int lps1(String str) {
		return lps1(str, 0, str.length() - 1);
	}

	public int lps1(String str, int i, int j) {
		if (i == j) return 1;
		if (str.charAt(i) == str.charAt(j) && i + 1 == j) return 2;
		if (str.charAt(i) == str.charAt(j)) return lps1(str, i + 1, j - 1) + 2;
		return Math.max(lps1(str, i, j - 1), lps1(str, i + 1, j));
	}

	// 3.DP-Bottom Up Approach
	public int lps3(String str) {
		int n = str.length();
		int[][] dp = new int[n][n];
		for (int i = 1; i < n; i++)
			dp[i][i] = 1;
		for (int len = 2; len <= n; len++) {
			for (int i = 0; i <= n - len; i++) {
				int j = i + len - 1;
				if (str.charAt(i) == str.charAt(j)) {
					dp[i][j] = len == 2 ? 2 : dp[i + 1][j - 1] + 2;
				} else {
					dp[i][j] = Math.max(dp[i][j - 1], dp[i + 1][j]);
				}
			}
		}

		printLPS3(dp, str);

		return dp[0][n - 1];
	}

	private String printLPS3(int[][] dp, String str) {
		int n = str.length(), maxLen = dp[0][n - 1];
		int row = 0, col = n - 1, i = 0, j = maxLen - 1;
		char[] seq = new char[maxLen];
		while (row <= col) {
			if (dp[row][col] > dp[row][col - 1] && dp[row][col] > dp[row + 1][col]) {
				seq[i++] = str.charAt(col);
				seq[j--] = str.charAt(col);
				row++;
				col--;
			} else if (dp[row][col] == dp[row][col - 1]) {
				col--;
			} else if (dp[row][col] == dp[row + 1][col]) {
				row++;
			} else {
				row++;
				col--;
			}
		}
		return String.valueOf(seq);
	}

	// Longest Palindromic Substring:
	/* Method 1(Brute	Force):
	 * The simple approach is to check each substring whether the substring is a palindrome or not. We can run three
	 * loops, the outer two loops pick all substrings one by one by fixing the corner characters, the inner loop checks
	 * whether the picked substring is palindrome or not.
	 * Time complexity: O ( n^3 )
	 */
	public String lpSubstr1(String str) {
		int max = -1, n = str.length();
		String maxString = null, subString;
		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				subString = str.substring(i, j + 1);
				if (isPalindrome(subString) && max < (j - i + 1)) {
					max = j - i + 1;
					maxString = subString;
				}
			}
		}
		return maxString;
	}

	// 2. Using DP-Bottom Up Approach:
	public String lpSubstr3(String str) {
		int n = str.length(), max = 1, start = 0;
		boolean[][] table = new boolean[n][n];
		for (int i = 0; i < n; i++)
			table[i][i] = true;
		for (int len = 2; len <= n; len++) {
			for (int i = 0; i <= n - len; i++) {
				int j = i + len - 1;
				if (str.charAt(i) == str.charAt(j) && (len == 2 || table[i + 1][j - 1])) {
					table[i][j] = true;
					if (len > max) {
						start = i;
						max = len;
					}
				}
			}
		}
		return str.substring(start, start + max);
	}

	public boolean isPalindrome(String str) {
		int l = 0, h = str.length() - 1;
		while (l < h) {
			if (str.charAt(l++) != str.charAt(h--)) return false;
		}
		return true;
	}

	// Palindromic Partitioning I
	public void palindromicPartitioningI(String s) {
		backtrackingPatterns.partition1(s);
		backtrackingPatterns.partition2(s);
	}

	/*
	 * Palindrome Partitioning II:
	 *   Given a string s, partition s such that every substring of the partition is a palindrome. Return the minimum cuts needed for 
	 *   a palindrome partitioning of s.
	 *   
	 *   Example: Input: s = "aab" Output: 1; 
	 *   Explanation: The palindrome partitioning ["aa","b"] could be produced using 1 cut.
	 */
	//Approach 1: Time O(n^2), Space: O(n^2)
	//Note: This is similar to palindromic substring problem.
	public int minCut1(String s) {
		int n = s.length();
		if (n <= 1) return 0;
		boolean[][] dp = new boolean[n][n];
		int[] cut = new int[n];

		for (int r = 0; r < n; r++) {
			int min = r;
			for (int l = 0; l <= r; l++) {
				if (s.charAt(l) == s.charAt(r) && (r - l <= 1 || dp[l + 1][r - 1])) {
					dp[l][r] = true;
					//if l == 0 means substring(0, r) is palindrome, so no cut is needed
					min = l == 0 ? 0 : Math.min(min, cut[l - 1] + 1);
				}
			}
			cut[r] = min;
		}
		return cut[n - 1];
	}

	//Approach 2: Similiar to above approach, except boolean array to mainatain palindrome
	//Time O(n^2), Space: O(n)
	public int minCut2(String s) {
		int n = s.length();
		int[] cuts = new int[n];

		for (int i = 0; i < n; i++)
			cuts[i] = i;

		for (int i = 0; i < n; i++) {
			checkPalindrome(s, cuts, i, i); //To handle odd length palindrome
			checkPalindrome(s, cuts, i, i + 1); //To handle even length palindrome
		}

		return cuts[n - 1];
	}

	private void checkPalindrome(String s, int[] cuts, int l, int r) {
		int n = cuts.length;
		while (l >= 0 && r < n && s.charAt(l) == s.charAt(r)) {
			//if l == 0 means substring(0, r) is palindrome, so no cut is needed
			if (l == 0) cuts[r] = 0;
			else cuts[r] = Math.min(cuts[r], cuts[l - 1] + 1);

			l--;
			r++;
		}
	}

	/************************** Pattern 6: String-Substring/Subsequence Probs *******************/
	// Longest Common Substring:
	// 1.Recursion Approach:
	public int lcStr1(String s1, String s2) {
		return lcStr1(s1, s2, s1.length() - 1, s2.length() - 1, 0);
	}

	public int lcStr1(String s1, String s2, int i, int j, int count) {
		if (i < 0 || j < 0) return count;
		if (s1.charAt(i) == s2.charAt(j)) return lcStr1(s1, s2, i - 1, j - 1, count + 1);
		return Utils.max(count, lcStr1(s1, s2, i - 1, j, 0), lcStr1(s1, s2, i, j - 1, 0));
	}

	// 2.DP:Bottom Up Approach:Time Complexity-O(m.n)
	public int lcStr3(String s1, String s2) {
		int m = s1.length(), n = s2.length();
		int[][] dp = new int[m][n];
		int max = 0, row = 0, col = 0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (s1.charAt(i) == s2.charAt(j)) {
					dp[i][j] = (i == 0 || j == 0) ? 1 : 1 + dp[i - 1][j - 1];
					if (max < dp[i][j]) {
						max = dp[i][j];
						row = i;
						col = j;
					}
				}
			}
		}
		//s1 - row or s2 - col
		printLCStr(dp, row, col, s1);
		return max;
	}

	public String printLCStr(int[][] dp, int row, int col, String s) {
		String subStr = "";
		while (row >= 0 && col >= 0 && dp[row][col] != 0) {
			subStr = s.charAt(row) + subStr;
			row--;
			col--;
		}
		return subStr;
	}

	// Longest Common subsequence:
	// 1.Recursive approach
	public int lcs1(String s1, String s2) {
		return lcs1(s1, s2, s1.length() - 1, s2.length() - 1);
	}

	private int lcs1(String s1, String s2, int i, int j) {
		if (i < 0 || j < 0) return 0;
		if (s1.charAt(i) == s2.charAt(j)) return 1 + lcs1(s1, s2, i - 1, j - 1);
		return Math.max(lcs1(s1, s2, i - 1, j), lcs1(s1, s2, i, j - 1));
	}

	// 2.DP:Top Down approach
	public int lcs2(String s1, String s2) {
		int[][] dp = new int[s1.length()][s2.length()];
		for (int[] row : dp)
			Arrays.fill(row, -1);
		return lcs2(s1, s2, s1.length() - 1, s2.length() - 1, dp);
	}

	private int lcs2(String s1, String s2, int i, int j, int[][] dp) {
		if (i < 0 || j < 0) return 0;
		if (dp[i][j] != -1) return dp[i][j];
		if (s1.charAt(i) == s2.charAt(j)) {
			return dp[i][j] = lcs2(s1, s2, i - 1, j - 1, dp) + 1;
		}
		return dp[i][j] = Math.max(lcs2(s1, s2, i - 1, j, dp), lcs2(s1, s2, i, j - 1, dp));
	}

	// 3.DP Bottom Up Approach
	public int lcs3(String s1, String s2) {
		int m = s1.length(), n = s2.length();
		int[][] dp = new int[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (s1.charAt(i) == s2.charAt(j)) {
					dp[i][j] = (i == 0 || j == 0) ? 1 : 1 + dp[i - 1][j - 1];
					dp[i][j] = 1 + dp[i - 1][j - 1];
				} else {
					dp[i][j] = Math.max(dp[i - 1][j], dp[i][j - 1]);
				}
			}
		}
		printLCS(dp, s1, s2);
		return dp[m - 1][n - 1];
	}

	// Print the longest common sub sequence
	private void printLCS(int[][] dp, String s1, String s2) {
		int i = s1.length(), j = s2.length();
		int longSeqCount = dp[i][j];
		char[] result = new char[longSeqCount];
		int index = longSeqCount;
		while (i >= 0 && j >= 0) {
			if (s1.charAt(i) == s2.charAt(j)) {
				result[--index] = s1.charAt(i);
				i--;
				j--;
			} else if (dp[i - 1][j] > dp[i][j - 1]) {
				i--;
			} else {
				j--;
			}
		}
		System.out.print("SubSequence:");
		for (int k = 0; k < longSeqCount; k++) {
			System.out.print(result[k] + "-");
		}
	}

	/* Edit Distance: Find minimum number of edits (operations) required to convert ‘str1’ into ‘str2’.
	 * You have the following three operations permitted on a word: 
	 * 	- Insert a character
	 * 	- Delete a character
	 * 	- Replace a character
	 * Input: word1 = "horse", word2 = "ros"; Output: 3
	 * Input: word1 = "intention", word2 = "execution"; Output: 5
	 */
	// Recursion Approach: O(3^n), n is Max(s1.length(), s2.length()).
	public int minDistance1(String s1, String s2) {
		return minDistance(s1, s2, s1.length() - 1, s2.length() - 1);
	}

	public int minDistance(String s1, String s2, int i, int j) {
		if (i < 0) return j + 1;
		if (j < 0) return i + 1;
		if (s1.charAt(i) == s2.charAt(j)) return minDistance(s1, s2, i - 1, j - 1);

		return 1 + Utils.min(minDistance(s1, s2, i, j - 1), //  represents insert operation
				minDistance(s1, s2, i - 1, j), // represents delete operation
				minDistance(s1, s2, i - 1, j - 1)); // represents replace operation
	}

	//DP- Top down Approach/Memoization: Time: O(mn), Space:O(mn)
	public int minDistance2(String s1, String s2) {
		int m = s1.length(), n = s2.length();
		int[][] memo = new int[m][n];
		for (int i = 0; i < m; i++)
			Arrays.fill(memo[i], -1);
		return minDistance2(s1, s2, m - 1, n - 1, memo);
	}

	public int minDistance2(String s1, String s2, int i, int j, int[][] memo) {
		if (i < 0) return j + 1;
		if (j < 0) return i + 1;
		if (memo[i][j] != -1) return memo[i][j];

		if (s1.charAt(i) == s2.charAt(j)) return memo[i][j] = minDistance2(s1, s2, i - 1, j - 1, memo);

		return memo[i][j] = 1 + Utils.min(minDistance2(s1, s2, i, j - 1, memo), minDistance2(s1, s2, i - 1, j, memo),
				minDistance2(s1, s2, i - 1, j - 1, memo));
	}

	// DP-Bottom up Approach: Time: O(mn), Space:O(mn)
	public int minDistance3(String s1, String s2) {
		int m = s1.length(), n = s2.length();
		if (m == 0 && n == 0) return 0;
		int[][] dp = new int[m + 1][n + 1];
		for (int i = 0; i <= m; i++) {
			for (int j = 0; j <= n; j++) {
				if (i == 0) dp[i][j] = j;
				else if (j == 0) dp[i][j] = i;
				else if (s1.charAt(i - 1) == s2.charAt(j - 1)) dp[i][j] = dp[i - 1][j - 1];
				else dp[i][j] = 1 + Utils.min(dp[i - 1][j - 1], dp[i - 1][j], dp[i][j - 1]);
			}
		}
		return dp[m][n];
	}

	/*
	 * Interleaving String:
	 *  Given s1, s2, s3, find whether s3 is formed by the interleaving of s1 and s2.
	 *  Example: 
	 *  	Input: s1 = "aabcc", s2 = "dbbca", s3 = "aadbbcbcac" Output: true
	 *  	Input: s1 = "aabcc", s2 = "dbbca", s3 = "aadbbbaccc" Output: false
	 */
	// Recursive Approach
	public boolean isInterleave1(String s1, String s2, String s3) {
		if ((s1.length() + s2.length()) != s3.length()) return false;
		return isInterleave(s1, s2, s3, 0, 0);
	}

	public boolean isInterleave(String s1, String s2, String s3, int i, int j) {
		if (i == s1.length() && j == s2.length() && i + j == s3.length()) return true;
		if (i + j == s3.length()) return false;

		return ((i < s1.length() && s1.charAt(i) == s3.charAt(i + j) && isInterleave(s1, s2, s3, i + 1, j))
				|| (j < s2.length() && s2.charAt(j) == s3.charAt(i + j) && isInterleave(s1, s2, s3, i, j + 1)));
	}

	// DP-Bottom Up Approach
	public boolean isInterleave3(String s1, String s2, String s3) {
		int n1 = s1.length(), n2 = s2.length();
		if ((n1 + n2) != s3.length()) return false;

		boolean[][] dp = new boolean[n1 + 1][n2 + 1];
		for (int i = 0; i <= n1; i++) {
			for (int j = 0; j <= n2; j++) {
				if (i == 0 && j == 0) {
					dp[i][j] = true;
				} else if (i == 0) {
					dp[i][j] = dp[i][j - 1] && s2.charAt(j - 1) == s3.charAt(i + j - 1);
				} else if (j == 0) {
					dp[i][j] = dp[i - 1][j] && s1.charAt(i - 1) == s3.charAt(i + j - 1);
				} else {
					dp[i][j] = (dp[i][j - 1] && s2.charAt(j - 1) == s3.charAt(i + j - 1))
							|| (dp[i - 1][j] && s1.charAt(i - 1) == s3.charAt(i + j - 1));
				}
			}
		}

		return dp[n1][n2];
	}

	/*
	 * Wild card Matching:
	 * Given an input string (s) and a pattern (p), implement wildcard pattern matching with support for '?' and '*'.
	'?' Matches any single character.
	'*' Matches any sequence of characters (including the empty sequence).
	 */
	// Approach1: Recursion
	public boolean wildCardMatch1(String s, String p) {
		return isMatch(0, s, 0, p);
	}

	private boolean isMatch(int i, String s, int j, String p) {
		int sn = s.length(), pn = p.length();
		if (j == pn && i == sn) return true;
		if (j == pn) return false;

		char pj = p.charAt(j);
		if (i < sn && (pj == '?' || pj == s.charAt(i))) {
			return isMatch(i + 1, s, j + 1, p);
		} else if (pj == '*') {
			return isMatch(i, s, j + 1, p) || i < sn && isMatch(i + 1, s, j, p);
		}

		return false;
	}

	// Approach2: Using DP-Top Down Approach
	public boolean wildCardMatch2(String s, String p) {
		Boolean[][] mem = new Boolean[s.length() + 1][p.length() + 1];
		return isMatch(0, s, 0, p, mem);
	}

	private boolean isMatch(int i, String s, int j, String p, Boolean[][] mem) {
		int sn = s.length(), pn = p.length();
		if (j == pn && i == sn) return true;
		if (j == pn) return false;

		if (mem[i][j] != null) return mem[i][j];

		char pj = p.charAt(j);
		if (i < sn && (pj == '?' || pj == s.charAt(i))) {
			return mem[i][j] = isMatch(i + 1, s, j + 1, p, mem);
		} else if (pj == '*') {
			return mem[i][j] = isMatch(i, s, j + 1, p, mem) || i < sn && isMatch(i + 1, s, j, p, mem);
		}

		return mem[i][j] = false;
	}

	// Approach3: Using DP-Bottom Up Approach- Time: O(mn), Space: O(mn)
	public boolean wildCardMatch3(String s, String p) {
		int m = s.length(), n = p.length();
		boolean[][] dp = new boolean[m + 1][n + 1];
		dp[0][0] = true;
		for (int j = 1; j <= n; j++)
			if (p.charAt(j - 1) == '*') dp[0][j] = dp[0][j - 1];

		for (int i = 1; i <= m; i++) {
			for (int j = 1; j <= n; j++) {
				//Here i-1, j-1 are curr index
				char si = s.charAt(i - 1), pj = p.charAt(j - 1);

				if (pj == si || pj == '?') {
					dp[i][j] = dp[i - 1][j - 1];
				} else if (pj == '*') {
					dp[i][j] = dp[i - 1][j] || dp[i][j - 1];
				} else {
					dp[i][j] = false;
				}

			}
		}
		return dp[m][n];
	}

	// Approach4: Greedy Algorithm
	/*
	 * Avg Time: O(m+n)
	 * Worst case definitely not linear. Should be O(mn)
	 * think that s ="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" to match p ="*aaaaaab"( '*' in the beginning)
	 * It's easy to see 'match' is moving step by step to almost the end, each time we move 'match', we will go 
	 * through the whole tail of p (after '*') until we found out 'b' is not a match. Thus it's O(NM)
	 */
	public int isMatch(final String s, final String p) {
		int m = s.length(), n = p.length();
		int i = 0, j = 0, star = -1, mark = 0;
		while (i < m) {
			if (j < n && (p.charAt(j) == s.charAt(i) || p.charAt(j) == '?')) {
				i++;
				j++;
			} else if (j < n && p.charAt(j) == '*') {
				star = j;
				mark = i;
				j++;
			} else if (star != -1) {
				j = star + 1;
				i = ++mark;
			} else {
				return 0;
			}
		}

		while (j < n && p.charAt(j) == '*') j++;

		return j == n ? 1 : 0;
	}

	/*
	 * Regular Expression Matching:
	 * Given an input string (s) and a pattern (p), implement regular expression matching with support for '.' and '*'.
	 *  '.' Matches any single character.
	 *  '*' Matches zero or more of the preceding element.
	 */
	// Approach1: Recursion
	public boolean regExMatch1(String s, String p) {
		return regEx(0, s, 0, p);
	}

	private boolean regEx(int i, String s, int j, String p) {
		int sn = s.length(), pn = p.length();
		if (j == pn && i == sn) return true;
		if (j == pn) return false;

		char pj = p.charAt(j);
		//Here first check "*", needs to look at the next char to repeat current char
		if (j + 1 < pn && p.charAt(j + 1) == '*') {
			if (regEx(i, s, j + 2, p)) return true;

			if (i < sn && (pj == '.' || pj == s.charAt(i))) {
				if (regEx(i + 1, s, j, p)) {
					return true;
				}
			}
		} else if (i < sn && (s.charAt(i) == pj || pj == '.')) {
			return regEx(i + 1, s, j + 1, p);
		}
		return false;
	}

	// Approach2: Using DP-Top Down Approach
	public boolean regExMatch2(String s, String p) {
		Boolean[][] mem = new Boolean[s.length() + 1][p.length() + 1];
		return regEx(0, s, 0, p, mem);
	}

	private boolean regEx(int i, String s, int j, String p, Boolean[][] mem) {
		int sn = s.length(), pn = p.length();
		if (j == pn && i == sn) return true;
		if (j == pn) return false;

		if (mem[i][j] != null) return mem[i][j];

		char pj = p.charAt(j);
		//Here first check "*", needs to look at the next char to repeat current char
		if (j + 1 < pn && p.charAt(j + 1) == '*') {
			if (regEx(i, s, j + 2, p, mem)) {
				return mem[i][j] = true;
			}
			if (i < sn && (pj == '.' || pj == s.charAt(i))) {
				if (regEx(i + 1, s, j, p, mem)) {
					return mem[i][j] = true;
				}
			}
		} else if (i < sn && (s.charAt(i) == pj || pj == '.')) {
			return mem[i][j] = regEx(i + 1, s, j + 1, p, mem);
		}
		return mem[i][j] = false;
	}

	// Approach3: Using DP-Bottom Up Approach- Time: O(mn), Space: O(mn) 
	public boolean regExMatch3(String s, String p) {
		int m = s.length(), n = p.length();
		boolean[][] dp = new boolean[m + 1][n + 1];
		// Base case: For both s & p are empty
		dp[0][0] = true;
		for (int j = 2; j <= n; j++)
			if (p.charAt(j - 1) == '*') dp[0][j] = dp[0][j - 2];

		for (int i = 1; i <= m; i++) {
			for (int j = 1; j <= n; j++) {
				// Note: Here i-1, j-1 is curr index/char;
				char si = s.charAt(i - 1), pj = p.charAt(j - 1);

				if (pj == si || pj == '.') {
					dp[i][j] = dp[i - 1][j - 1];
				} else if (pj == '*') {
					dp[i][j] = dp[i][j - 2]; //This is for 0 occurrence of the prev char
					//Check curr char in 's' with previous char in pattern 'p'
					if (p.charAt(j - 2) == si || p.charAt(j - 2) == '.') {
						dp[i][j] = dp[i][j] || dp[i - 1][j];
					}
				} else {
					dp[i][j] = false;
				}
			}
		}
		return dp[m][n];
	}
}