package com.basic.algorithms;

import com.basic.algorithms.operations.PatternSearchOperations;

/*
 * The Pattern Searching algorithms are sometimes also referred to as String Searching Algorithms and are considered 
 * as a part of the String algorithms. These algorithms are useful in the case of searching a string within another string.
 * 
 * Pattern searching is an important problem in computer science. When we do search for a string in notepad/word
 * file or browser or database, pattern searching algorithms are used to show the search results.
 */
public class PatternSearchingAlgorithms implements PatternSearchOperations {

	/* Naive Pattern Searching:- Time complexity:O(mn) = Exactly O(m(n-m))
	 * The Naive String Matching algorithm slides the pattern one by one. After each slide, it one by one checks characters at the current
	 * shift and if all characters match then prints the match.
	 * Best case: O(n)
	 * 	Eg: txt = "AABCCAADDEE", pat = "FAA";
	 * Worst Case: O(n*(n-m));
	 * 	Eg: txt = "AAAAAAAAAAAAAAAAAA"; pat = "AAAAA"; 
	 */
	@Override
	public void naivePatternSearching(String txt, String pattern) {
		int j;
		int n = txt.length(); // String Length
		int m = pattern.length(); // Pattern Length
		for (int i = 0; i <= (n - m); i++) {
			if (txt.charAt(i) != pattern.charAt(0)) continue;
			for (j = 0; j < m; j++) {
				if (txt.charAt(i + j) != pattern.charAt(j)) break;
			}
			if (j == m) System.out.println("Pattern found at: " + i);
		}
	}

	/* KMP(Knuth Morris Pratt) Algorithm; Time=O(n+m); Space: O(n)
	 * The KMP matching algorithm uses degenerating property (pattern having same sub-patterns appearing more than once in the pattern) 
	 * of the pattern and improves the worst case complexity to O(n). The basic idea behind KMP’s algorithm is: whenever we detect a 
	 * mismatch (after some matches), we already know some of the characters in the text of the next window. We take advantage of this 
	 * information to avoid matching the characters that we know will anyway match. 
	 */
	@Override
	public void KMPAlgorithm(String text, String pattern) {
		int m = pattern.length(), n = text.length();
		int[] lps = computeLPSArray(pattern);

		//Here 'i' is txt index, j is pat index 
		int i = 0, j = 0;

		while (i < n) {
			if (text.charAt(i) == pattern.charAt(j)) {
				i++;
				j++;

				if (j == m) {
					System.out.println("Pattern found at index:" + (i - j));
					//move to next highest prefix in the pattern
					j = lps[j - 1];
				}
			} else if (j != 0) {
				//Here if char is not match and index j is not zero, then move to next highest prefix in the pattern 
				j = lps[j - 1];
			} else {
				i++;
			}
		}
	}

	/* 
	 * Preprocessing Overview:
	 * - KMP algorithm preprocesses pat[] and constructs an auxiliary lps[] of size m (same as size of pattern) which is used to
	 *   skip characters while matching.
	 * - Name lps indicates longest proper prefix which is also suffix.. A proper prefix is prefix with whole string not allowed.
	 *   For example, prefixes of “ABC” are “”, “A”, “AB” and “ABC”. Proper prefixes are “”, “A” and “AB”. Suffixes of the string 
	 *   are “”, “C”, “BC” and “ABC”.
	 * - We search for lps in sub-patterns. More clearly we focus on sub-strings of patterns that are either prefix and suffix.
	 * - For each sub-pattern pat[0..i] where i = 0 to m-1, lps[i] stores length of the maximum matching proper prefix which is 
	 *   also a suffix of the sub-pattern pat[0..i]. 
	 *   lps[i] = the longest proper prefix of pat[0..i] which is also a suffix of pat[0..i]. 
	 *   
	 *   
	 * How to use lps[] to decide next positions (or to know a number of characters to be skipped)? 
	 *  1.We start comparison of pat[j] with j = 0 with characters of current window of text.
	 *  2.We keep matching characters txt[i] and pat[j] and keep incrementing i and j while pat[j] and txt[i] keep matching.
	 *  3.When we see a mismatch
	 *      - We know that characters pat[0..j-1] match with txt[i-j…i-1] (Note that j starts with 0 and increment it only when there is a match).
	 *      - We also know (from above definition) that lps[j-1] is count of characters of pat[0…j-1] that are both proper prefix and suffix.
	 *      - From above two points, we can conclude that we do not need to match these lps[j-1] characters with txt[i-j…i-1] because we know that 
	 *        these characters will anyway match. Let us consider above example to understand this.
	 */
	public int[] computeLPSArray(String pat) {
		int m = pat.length();
		int[] lps = new int[m]; //LPS - Longest Prefix Suffix:
		lps[0] = 0;
		//Here we use two ptrs l & r to preprocess the pattern string
		int l = 0, r = 1;
		while (r < m) {
			if (pat.charAt(r) == pat.charAt(l)) {
				// l= increment & assign, because we can start checking from the next index in pattern, when we search pattern in text
				lps[r++] = ++l;
			} else if (l != 0) {
				/* If index l is not zero, then we will move to next highest prefix, to check whether same prefix and suffix 
				 * available in the pattern */
				l = lps[l - 1];
			} else {
				lps[r] = 0;
				r++;
			}

		}
		return lps;
	}

	/* Rabin-Karp Algorithm: (Rolling Hash Technique)
	 * The Naive String Matching algorithm slides the pattern one by one. After each slide, it one by one checks characters at the current
	 * shift and if all characters match then prints the match. 
	 * Like the Naive Algorithm, Rabin-Karp algorithm also slides the pattern one by one. But unlike the Naive algorithm, Rabin Karp algorithm
	 * matches the hash value of the pattern with the hash value of current substring of text, and if the hash values match then only it starts
	 * matching individual characters. 
	 * Note: abin Karp Algorithm works based on Rolling Hash Algorithm. So we can use any logic to calculate the hash values for following strings.
	 * 	1) Pattern itself (Hashing)
	 * 	2) All the substrings of the text of length m (Hashing & Rehashing)
	 * 
	 * The average and best-case running time of the Rabin-Karp algorithm is O(n+m), but its worst-case time is O(nm).
	 * 
	 */
	@Override
	public void rabinKarpAlgorithm(String text, String pat) {
		int prime = 101; //It can be 11, 101, etc
		long textHash = 0, patternHash = 0;
		int n = text.length(), m = pat.length();

		// Calculate hash value of pat and text: index 0 to m-1
		for (int i = 0; i < m; i++) {
			// Pattern Hash Value
			patternHash += pat.charAt(i) * Math.pow(prime, i);
			// Hash Value of first m character in Text
			textHash += text.charAt(i) * Math.pow(prime, i);
		}
		int j;
		for (int i = 0; i <= n - m; i++) {
			if (textHash == patternHash) {
				// Check for characters one by one
				for (j = 0; j < m; j++) {
					if (text.charAt(i + j) != pat.charAt(j)) break;
				}
				if (j == m) System.out.println("Pattern found at index:" + i);
			}

			// ReHashing
			// Note: Prime value's index always from index 0 to m-1;
			if (i < n - m) {
				textHash -= text.charAt(i); // X = OldHash - text(first char); 
				textHash /= prime; // X = X/prime
				textHash += (text.charAt(i + m) * Math.pow(prime, m - 1)); // NewHash = X + text(next char) *
																			// prime^(m-1);
			}
		}
	}

	/*
	 * Approach2:
	 *   Rabin Karp algorithm is same as above, to check the hash value of the pattern and text. 
	 *   Here we use different logic to create the hashing and rehashing. This logic avoids integer overflow.
	 *   
	 *   TODO: Revisit this approach
	 *   Reference: 
	 *   	https://www.youtube.com/watch?v=d3TZpfnpJZ0
	 *   	https://www.youtube.com/watch?v=BQ9E-2umSWc
	 *   	https://www.geeksforgeeks.org/rabin-karp-algorithm-for-pattern-searching/
	 */

	/* d is the number of characters in the input alphabet. It may be 26(Only lower case alphabets), 52(both lower and upper case alphabets),
	 * 128(ASCII Chars) or 256(All ASCII Chars)*/
	public final static int d = 256;

	public void rabinKarpAlgorithm2(String txt, String pat) {
		int M = pat.length();
		int N = txt.length();
		int prime = 101; // Prime Number
		int i, j;
		int patHash = 0; // hash value for pattern
		int txtHash = 0; // hash value for txt
		int h = 1;

		// The value of h would be "pow(d, M-1)%prime"
		for (i = 0; i < M - 1; i++)
			h = (h * d) % prime;

		// Calculate the hash value of pattern and first
		// window of text
		for (i = 0; i < M; i++) {
			patHash = (d * patHash + pat.charAt(i)) % prime;
			txtHash = (d * txtHash + txt.charAt(i)) % prime;
		}

		// Slide the pattern over text one by one
		for (i = 0; i <= N - M; i++) {

			// Check the hash values of current window of text
			// and pattern. If the hash values match then only
			// check for characters on by one
			if (patHash == txtHash) {
				/* Check for characters one by one */
				for (j = 0; j < M; j++) {
					if (txt.charAt(i + j) != pat.charAt(j)) break;
				}

				// if patHash == txtHash and pat[0...M-1] = txt[i, i+1, ...i+M-1]
				if (j == M) System.out.println("Pattern found at index " + i);
			}

			// Calculate hash value for next window of text: Remove
			// leading digit, add trailing digit
			if (i < N - M) {
				txtHash = (d * (txtHash - txt.charAt(i) * h) + txt.charAt(i + M)) % prime;

				// We might get negative value of txtHash, converting it
				// to positive
				if (txtHash < 0) txtHash = (txtHash + prime);
			}
		}
	}

	@Override
	public void finiteAutomata(String str, String pattern) {
		// TODO Auto-generated method stub

	}

	// Z Algorithm
	@Override
	public void zAlgorithm(String text, String pattern) {
		int zArrayLength = text.length() + pattern.length() + 1;
		String zString = pattern + "$" + text;
		int[] z = calculateZArray(zString);

		for (int i = 0; i < zArrayLength; i++) {
			if (z[i] == pattern.length()) {
				System.out.println("Pattern found at index: " + (i - pattern.length() - 1));
			}
		}
	}

	private int[] calculateZArray(String str) {
		int n = str.length();
		int[] zArray = new int[n];
		int left = 0, right = 0;
		int j;
		for (int i = 1; i < n; i++) {
			if (i > right) {
				left = right = i;
				while (right < n && str.charAt(right - left) == str.charAt(right)) right++;

				zArray[i] = right - left;
				right--;
			} else {
				// j corresponds to number which matches in [L,R] interval.
				j = i - left;

				/*if Z[j] is less than remaining interval then Z[i] will be equal to Z[j]. For example, str = "ababab", i = 3, R = 5 and L = 2*/
				if (zArray[j] < right - 1 + 1) {
					zArray[i] = zArray[j];
				} else {
					// else start from R and check manually
					left = i;
					while (right < n && str.charAt(i) == str.charAt(right - left)) right++;

					zArray[i] = right - left;
					right--;
				}
			}
		}
		return zArray;
	}

}