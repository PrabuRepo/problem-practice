package com.problems.patterns.crossdomains;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.common.model.TrieNode;
import com.common.model.WordNode;

public class WordProblems {

	KElementsPattern kElementsPattern;

	/* Shortest Word Distance:
	 * words = ["practice", "makes", "perfect", "coding", "makes"]. 
	 * Given word1 = "coding", word2 = "practice", Result: 3. 
	 * Given word1 = "makes", word2 = "coding", Result: 1.
	 */
	public int shortestDistanceI(String[] words, String word1, String word2) {
		int min = Integer.MAX_VALUE, index1 = -1, index2 = -1;
		for (int i = 0; i < words.length; i++) {
			if (words[i].equals(word1)) index1 = i;

			if (words[i].equals(word2)) index2 = i;

			if (index1 != -1 && index2 != -1) min = Math.min(min, Math.abs(index1 - index2));
		}
		return min;
	}

	/* Word Pattern: Example 1: Input: pattern = "abba", 
	 * str = "dog cat cat dog";  Output: true
	 */
	//Approach1: Using Map & Set - Time: O(n), Space: O(2n) = O(n)
	public boolean wordPattern1(String pattern, String str) {
		if (str == null || str.length() == 0 || pattern.length() == 0) return false;

		Map<Character, String> map = new HashMap<>();
		Set<String> set = new HashSet<>();
		String[] words = str.split(" ");

		if (pattern.length() != words.length) return false;

		for (int i = 0; i < words.length; i++) {
			char ch = pattern.charAt(i);
			if (map.containsKey(ch)) {
				if (!map.get(ch).equals(words[i])) return false;
			} else {
				//Validation to check any duplicate word present in the words array
				if (!set.add(words[i])) return false;

				map.put(ch, words[i]);
			}
		}

		return true;
	}

	//Approach2: Using two Maps to store the index - Time: O(n), Space: O(2n) = O(n)
	public boolean wordPattern2(String pattern, String str) {
		if (str == null || str.length() == 0 || pattern.length() == 0) return false;

		Map<Character, Integer> charsMap = new HashMap<>();
		Map<String, Integer> wordsMap = new HashMap<>();
		String[] words = str.split(" ");

		if (pattern.length() != words.length) return false;

		for (int i = 0; i < words.length; i++) {
			char ch = pattern.charAt(i);
			if (!charsMap.containsKey(ch)) charsMap.put(ch, i);
			if (!wordsMap.containsKey(words[i])) wordsMap.put(words[i], i);

			if (!charsMap.get(ch).equals(wordsMap.get(words[i]))) return false;
		}

		return true;
	}

	/*
	 * Word Pattern II: Given a pattern and a string str, find if str follows the same pattern.
	 * Here non-empty substring in str.
	 * 	Examples:	
	 * 		pattern = "cdcd", str = "xyabxyab" should return true.
	 * 		pattern = "abab", str = "redblueredblue" should return true.
	 * 		pattern = "aaaa", str = "asdasdasdasd" should return true.
	 * 		pattern = "aabb", str = "xyzabcxzyabc" should return false.
	 */
	public boolean wordPatternMatch1(String pattern, String str) {
		if (pattern.length() == 0 && str.length() == 0) return true;
		if (pattern.length() == 0) return false;

		HashMap<Character, String> map = new HashMap<>();
		Set<String> set = new HashSet<>();
		return helper(pattern, str, 0, 0, map, set);
	}

	/*
	 * Solution:
	 *   Iterate the chars in str and Apply the Word Break I solution, if substr doesn't matches pattern
	 *   remove substr from map and proceed with the next one. 
	 */
	public boolean helper(String pattern, String str, int i, int j, Map<Character, String> map, Set<String> set) {
		if (i == pattern.length() && j == str.length()) return true;

		if (i >= pattern.length() || j >= str.length()) return false;

		char pat = pattern.charAt(i); // pattern char
		for (int k = j + 1; k <= str.length(); k++) {
			String sub = str.substring(j, k);
			if (!map.containsKey(pat) && !set.contains(sub)) {
				map.put(pat, sub);
				set.add(sub);

				if (helper(pattern, str, i + 1, k, map, set)) return true;

				// Backtracking, remove and check from next index
				map.remove(pat);
				set.remove(sub);
			} else if (map.containsKey(pat) && map.get(pat).equals(sub)) {
				if (helper(pattern, str, i + 1, k, map, set)) return true;
			}
		}

		return false;
	}

	/* Since containsValue() method is used here, the time complexity is O(n). We can use another set to track the value 
	 * set which leads to time complexity of O(1):
	 */
	public boolean wordPatternMatch2(String pattern, String str) {
		if (pattern.length() == 0 && str.length() == 0) return true;
		if (pattern.length() == 0) return false;

		HashMap<Character, String> map = new HashMap<>();
		HashSet<String> set = new HashSet<>();
		return helper(pattern, str, 0, 0, map, set);
	}

	public boolean helper(String pattern, String str, int i, int j, HashMap<Character, String> map,
			HashSet<String> set) {
		if (i == pattern.length() && j == str.length()) return true;

		if (i >= pattern.length() || j >= str.length()) return false;

		char c = pattern.charAt(i);
		for (int k = j + 1; k <= str.length(); k++) {
			String sub = str.substring(j, k);
			if (!map.containsKey(c) && !set.contains(sub)) {
				map.put(c, sub);
				set.add(sub);
				if (helper(pattern, str, i + 1, k, map, set)) return true;

				// Backtracking, remove and check from next index
				map.remove(c);
				set.remove(sub);
			} else if (map.containsKey(c) && map.get(c).equals(sub)) {
				if (helper(pattern, str, i + 1, k, map, set)) return true;
			}
		}

		return false;
	}

	/*
	 * Word Break I:
	 * Given a non-empty string s and a dictionary wordDict containing a list of non-empty words, determine if s can
	 * be segmented into a space-separated sequence of one or more dictionary words.
	 * 	Note:The same word in the dictionary may be reused multiple times in the segmentation. You may assume the 
	 * dictionary does not contain duplicate words.
	 * 	Example 1:	Input: s = "leetcode", wordDict = ["leet", "code"]	Output: true
	 * 	Explanation: Return true because "leetcode" can be segmented as "leet code".
	 * 	Input: s = "catsandog", wordDict = ["cats", "dog", "sand", "and", "cat"] Output: false 
	 */

	/* Solution:
	 * This problem cane be solved by two ways,
	 *   1.By iterating words in the dictionary and check whether word is present in string 
	 *       - wordBreakI11 -> Bactracking solution
	 *       - wordBreakI21 -> DP solution
	 *   2.By iterating chars in the string and check whether substring is present in the dictionary.
	 *       - wordBreakI12 -> Backtracking solution
	 *       - wordBreakI22 -> DP solution 
	 */
	// Using DFS Search: It throws TLE
	public boolean wordBreakI11(String s, List<String> wordDict) {
		return wordBreakHelper(s, wordDict, 0);
	}

	private boolean wordBreakHelper(String s, List<String> dict, int index) {
		if (index == s.length()) return true;

		for (String word : dict) {
			int len = index + word.length();

			if (len > s.length()) continue;

			String substr = s.substring(index, len);
			if (substr.equals(word)) {
				if (wordBreakHelper(s, dict, len)) return true;
			}
		}

		return false;
	}

	public boolean wordBreakI12(String s, List<String> wordDict) {
		return wordBreakHelper(s, new HashSet<String>(wordDict), 0);
	}

	public boolean wordBreakHelper(String s, Set<String> dict, int index) {
		if (index == s.length()) return true;

		for (int i = index + 1; i <= s.length(); i++) {
			String sub = s.substring(index, i);
			if (!dict.contains(sub)) continue;

			if (wordBreakHelper(s, dict, i)) return true;
		}

		return false;
	}

	// DP: Using string length & dict size; Time: O(string length * dict size).
	public boolean wordBreakI21(String s, List<String> wordDict) {
		int n = s.length();
		boolean[] lookup = new boolean[n + 1];
		lookup[0] = true;

		for (int i = 0; i < n; i++) {
			if (!lookup[i]) continue;
			for (String word : wordDict) {
				int len = i + word.length();
				if (len > s.length() || lookup[len]) continue;
				if (s.substring(i, len).equals(word)) {
					lookup[len] = true;
				}
			}
		}
		return lookup[n];
	}

	// DP: Using only string length; Time: O(string length * string length); Space:O(n+dictSize)
	public boolean wordBreakI22(String s, List<String> wordDict) {
		int n = s.length();
		boolean[] lookup = new boolean[n + 1];
		lookup[0] = true;
		Set<String> set = new HashSet<>();
		set.addAll(wordDict);

		for (int i = 0; i < n; i++) {
			if (!lookup[i]) continue;
			for (int j = i + 1; j <= n; j++)
				if (set.contains(s.substring(i, j))) {
					lookup[j] = true;
				}
		}

		return lookup[n];
	}

	/*
	 * Word Break II:
	 * Return all such possible sentences.
	 * Example 1: 
	 *  Input: s = "pineapplepenapple"; wordDict = ["apple", "pen", "applepen", "pine", "pineapple"]
	 *  Output:["pine apple pen apple", "pineapple pen apple","pine applepen apple"]
	 *  Explanation: Note that you are allowed to reuse a dictionary word.
	 */
	// Using DFS Search: It throws TLE
	public List<String> wordBreakII1(String s, List<String> wordDict) {
		List<String> result = new ArrayList<>();
		wordBreakHelper(s, wordDict, 0, "", result);
		return result;
	}

	private void wordBreakHelper(String s, List<String> dict, int index, String str, List<String> result) {
		if (index == s.length()) {
			result.add(str.trim());
			return;
		}

		for (String word : dict) {
			int end = index + word.length();

			if (end > s.length()) continue;

			String substr = s.substring(index, end);
			if (substr.equals(word)) {
				wordBreakHelper(s, dict, end, str + " " + substr, result);
			}
		}
	}

	// DP: Using only string length; Time: O(string length * string length).
	public List<String> wordBreakII2(String s, List<String> wordDict) {
		int n = s.length();
		List<String>[] lookup = new ArrayList[n + 1];
		lookup[0] = new ArrayList<>();

		for (int i = 0; i < n; i++) {
			if (lookup[i] == null) continue;
			for (int j = i + 1; j <= n; j++) {
				String subStr = s.substring(i, j);
				if (wordDict.contains(subStr)) {
					if (lookup[j] == null) lookup[j] = new ArrayList<>();
					lookup[j].add(subStr);
				}
			}

		}

		List<String> result = new ArrayList<>();
		if (lookup[n] == null) return result;
		dfs(lookup, result, "", s.length());

		return result;
	}

	private void dfs(List<String>[] lookup, List<String> result, String str, int i) {
		if (i == 0) {
			result.add(str.trim());
		} else {
			for (String word : lookup[i])
				dfs(lookup, result, word + " " + str, i - word.length());
		}
	}

	/*
	 * Valid Word Square:
	 * Given a sequence of words, check whether it forms a valid word square.
	 * For example, the word sequence ["ball","area","lead","lady"] forms a word square because each word reads the same both horizontally and vertically.
	 * b a l l
	 * a r e a
	 * l e a d
	 * l a d y
	 */
	public boolean validWordSquare(List<String> words) {
		if (words == null || words.size() == 0) return true;
		int m = words.size();
		for (int i = 0; i < m; i++) {
			int n = words.get(i).length();
			for (int j = 0; j < n; j++)
				if (j >= m || m != n || words.get(i).charAt(j) != words.get(j).charAt(i)) return false;
		}

		return true;
	}

	/*
	 * Word Squares:
	 * Given a set of words (without duplicates), find all word squares you can build from them.
	 * A sequence of words forms a valid word square if the kth row and column read the exact same string, where 0 <=k < max(numRows, numColumns).
	 */
	TrieNode root;

	public List<List<String>> wordSquares(String[] words) {
		List<List<String>> ans = new ArrayList<>();
		if (words == null || words.length == 0) return ans;
		int len = words[0].length();
		// Build Trie
		buildTrie(words);
		// Search Words
		List<String> ansBuilder = new ArrayList<>();
		for (String w : words) {
			ansBuilder.add(w);
			search(len, ans, ansBuilder);
			ansBuilder.remove(ansBuilder.size() - 1);
		}
		ans.stream().forEach(k -> System.out.print(k + ", "));
		return ans;
	}

	void buildTrie(String[] words) {
		root = new TrieNode();
		for (String w : words) {
			TrieNode cur = root;
			for (char ch : w.toCharArray()) {
				int idx = ch - 'a';
				if (cur.children[idx] == null) cur.children[idx] = new TrieNode();
				cur.children[idx].startWith.add(w);
				cur = cur.children[idx];
			}
		}
	}

	private void search(int len, List<List<String>> ans, List<String> ansBuilder) {
		if (ansBuilder.size() == len) {
			ans.add(new ArrayList<>(ansBuilder));
			return;
		}

		int idx = ansBuilder.size();
		StringBuilder prefixBuilder = new StringBuilder();
		for (String s : ansBuilder)
			prefixBuilder.append(s.charAt(idx));
		List<String> startWith = findByPrefix(prefixBuilder.toString());
		for (String sw : startWith) {
			ansBuilder.add(sw);
			search(len, ans, ansBuilder);
			ansBuilder.remove(ansBuilder.size() - 1);
		}
	}

	List<String> findByPrefix(String prefix) {
		List<String> ans = new ArrayList<>();
		TrieNode cur = root;
		for (char ch : prefix.toCharArray()) {
			int idx = ch - 'a';
			if (cur.children[idx] == null) return ans;

			cur = cur.children[idx];
		}
		ans.addAll(cur.startWith);
		return ans;
	}

	private static final int[][] DIRS = { { 0, 1 }, { 0, -1 }, { 1, 0 }, { -1, 0 } };

	//Matrix 4 Dir & Trie: Word Boggle or Word Search I, II
	/* Word Search I - Search one word
	 * Given a 2D board and a "word", find if the word exists in the grid.The word can be constructed from 
	 * letters of sequentially adjacent cell, where "adjacent" cells are those horizontally or vertically neighboring. 
	 * The same letter cell may not be used more than once.
	 * Example:
	 * board =[['A','B','C','E'], ['S','F','C','S'], ['A','D','E','E']]
	 * Given word = "ABCCED", return true.
	 * Given word = "SEE", return true.
	 */
	/*
	 * Time is O(M * N * 4^L) where M*N is the size of the board and we have 4^L for each cell because of the recursion. 
	 * where L is the length of the word;
	 */
	public boolean wordSearchI(char[][] board, String str) {
		if (str.length() == 0 || board.length == 0 || board[0].length == 0) return false;

		int row = board.length, col = board[0].length;

		for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++)
				if (str.charAt(0) == board[i][j]) {
					if (dfsSearch(board, str, i, j, 0)) return true;
				}

		return false;
	}

	private boolean dfsSearch(char[][] board, String word, int i, int j, int index) {
		int row = board.length, col = board[0].length;
		if (i < 0 || i >= row || j < 0 || j >= col || index >= word.length() || word.charAt(index) != board[i][j]
				|| board[i][j] == '#')
			return false;
		if (index == word.length() - 1) return true;

		char temp = board[i][j];
		board[i][j] = '#'; // Avoid to revisit the same value

		boolean flag = false;
		flag = dfsSearch(board, word, i - 1, j, index + 1) || dfsSearch(board, word, i + 1, j, index + 1)
				|| dfsSearch(board, word, i, j - 1, index + 1) || dfsSearch(board, word, i, j + 1, index + 1);

		//or
		for (int[] dir : DIRS) {
			flag = flag || dfsSearch(board, word, i + dir[0], j + dir[1], index + 1);
		}

		board[i][j] = temp;
		return flag;
	}

	/* Word Search II - Search array of words: 
	 * Given a 2D board and a list of words from the dictionary, find all words in the board.
	 * Each word must be constructed from letters of sequentially adjacent cell, where "adjacent" cells are those horizontally 
	 * or vertically neighboring. The same letter cell may not be used more than once in a word.
	 * 	Example:
	 * 	Input: words = ["oath","pea","eat","rain"] and board =
	 * 	[ ['o','a','a','n'],
	 *    ['e','t','a','e'],
	 *    ['i','h','k','r'],
	 *    ['i','f','l','v']
	 *  ]
	 *  Output: ["eat","oath"]
	 */
	// Approach1: Using DFS -> Time Complexity: O(len*m*n) where len- no of words, m- row size, n-colSize
	public List<String> wordSearchII1(char[][] board, String[] words) {
		List<String> result = new ArrayList<>();
		if (words.length == 0 || board.length == 0 || board[0].length == 0) return result;

		HashSet<String> set = new HashSet<>(); // Set is used to remove the duplicate word
		for (String word : words)
			if (!set.add(word) && wordSearchI(board, word)) {
				result.add(word);
			}

		result.stream().forEach(k -> System.out.print(k + " "));
		return result;
	}

	public List<String> wordSearchII2(char[][] board, String[] words) {
		List<String> result = new ArrayList<>();

		// Build Trie data structure
		TrieNode root = new TrieNode();
		for (String word : words) {
			insert1(root, word);
		}

		// dfs search
		for (int i = 0; i < board.length; i++)
			for (int j = 0; j < board[0].length; j++)
				dfsSearch1(board, root, i, j, result);

		result.stream().forEach(k -> System.out.print(k + " "));
		return result;
	}

	public void dfsSearch1(char[][] board, TrieNode root, int i, int j, List<String> result) {
		int rSize = board.length, cSize = board[0].length;
		// Row & col Validation
		if (i < 0 || i >= rSize || j < 0 || j >= cSize) return;
		// Trie Validation
		char ch = board[i][j];
		if (ch == '#' || root.childNodes.get(ch) == null) return;

		root = root.childNodes.get(ch);
		//Record and move on
		if (root.word != null) {
			result.add(root.word);
			root.word = null;
		}

		board[i][j] = '#';

		dfsSearch1(board, root, i, j - 1, result);
		dfsSearch1(board, root, i, j + 1, result);
		dfsSearch1(board, root, i - 1, j, result);
		dfsSearch1(board, root, i + 1, j, result);

		board[i][j] = ch;
	}

	public List<String> wordSearchII3(char[][] board, String[] words) {
		List<String> result = new ArrayList<>();

		// Build Trie data structure
		TrieNode root = new TrieNode();
		for (String word : words) {
			insert2(root, word);
		}

		// dfs search
		for (int i = 0; i < board.length; i++)
			for (int j = 0; j < board[0].length; j++)
				dfsSearch2(board, root, i, j, new StringBuilder(), result);

		result.stream().forEach(k -> System.out.print(k + " "));
		return result;
	}

	public void dfsSearch2(char[][] board, TrieNode root, int i, int j, StringBuilder sb, List<String> result) {
		int rSize = board.length, cSize = board[0].length;
		// Row & col Validation
		if (i < 0 || i >= rSize || j < 0 || j >= cSize) return;
		// Trie Validation
		char ch = board[i][j];
		if (ch == '#' || root.childNodes.get(ch) == null) return;

		root = root.childNodes.get(ch);

		if (root.isEndOfWord) {
			result.add(sb.toString() + ch);
			root.isEndOfWord = false;
			return;
		}

		sb.append(ch);
		board[i][j] = '#';

		dfsSearch2(board, root, i, j - 1, sb, result);
		dfsSearch2(board, root, i, j + 1, sb, result);
		dfsSearch2(board, root, i - 1, j, sb, result);
		dfsSearch2(board, root, i + 1, j, sb, result);

		board[i][j] = ch;
	}

	// Insert all the words in the Trie DS
	public void insert1(TrieNode root, String word) {
		TrieNode curr = root;
		for (char ch : word.toCharArray()) {
			curr = curr.childNodes.computeIfAbsent(ch, node -> new TrieNode());
			/*
			if (curr.childNodes.get(ch) == null) curr.childNodes.put(ch, new TrieNode());
			curr = curr.childNodes.get(ch);*/
		}
		curr.word = word;
	}

	// Insert all the words in the Trie DS
	public void insert2(TrieNode root, String word) {
		TrieNode curr = root;
		for (char ch : word.toCharArray()) {
			curr = curr.childNodes.computeIfAbsent(ch, node -> new TrieNode());
		}
		curr.isEndOfWord = true;
	}

	/* Word Ladder: Find the ladder length
	 * Given two words (beginWord and endWord), and a dictionary's word list, find the length of shortest transformation sequence 
	 * from beginWord to endWord, such that:
	 *    Only one letter can be changed at a time.
	 *    Each transformed word must exist in the word list. Note that beginWord is not a transformed word.
	 *    Input: beginWord = "hit", endWord = "cog"; wordList = ["hot","dot","dog","lot","log","cog"]; 
	 *    Output: 5;
	 *    Explanation: As one shortest transformation is "hit" -> "hot" -> "dot" -> "dog" -> "cog", return its length 5.
	 */
	//Approach1: Using BFS and WordNode
	/*
	 * Time complexity: O(N * 26 * wordLength^2) time. where N=wordList.size(), M = Word length
	 * For each word in the word list, we iterate over its length to find all the intermediate words corresponding to it. Since the length 
	 * of each word is M and we have N words, the total number of iterations the algorithm takes M×N. Additionally, forming each of the 
	 * intermediate word takes O(M) time because of the substring operation used to create the new string. This adds up to a complexity of O(M^2*N).
	 * Space complexity: O(wordList.size()) space. Because we add all words into a HashSet, and queue and seen set can't have more than wordList.size() elements.
	 */
	public int wordLadderI1(String beginWord, String endWord, List<String> wordList) {
		LinkedList<WordNode> queue = new LinkedList<>();
		queue.add(new WordNode(beginWord, 1));
		Set<String> wordDict = new HashSet<>(wordList);

		while (!queue.isEmpty()) {
			WordNode curr = queue.poll();
			if (curr.word.equals(endWord)) return curr.count;

			char[] arr = curr.word.toCharArray();
			for (int i = 0; i < arr.length; i++) {
				char temp = arr[i];
				for (char ch = 'a'; ch <= 'z'; ch++) {
					if (arr[i] == ch) continue;
					arr[i] = ch;
					String newStr = new String(arr);
					if (wordDict.contains(newStr)) {
						queue.add(new WordNode(newStr, curr.count + 1));
						wordDict.remove(newStr);
					}
					arr[i] = temp;
				}
			}
		}
		return 0;
	}

	//Approach2: Using BFS and String
	public int wordLadderI2(String beginWord, String endWord, List<String> wordList) {
		if (beginWord == null || endWord == null || wordList == null || wordList.size() == 0) return 0;

		LinkedList<String> queue = new LinkedList<>();
		Set<String> wordDict = new HashSet<>(wordList);
		queue.add(beginWord);
		int ladderLen = 1;

		while (!queue.isEmpty()) {
			int size = queue.size();
			while (size-- > 0) {
				String top = queue.poll();
				if (top.equals(endWord)) return ladderLen;
				char[] arr = top.toCharArray();
				for (int i = 0; i < arr.length; i++) {
					for (char c = 'a'; c <= 'z'; c++) {
						if (arr[i] == c) continue;
						char temp = arr[i];
						arr[i] = c;
						String newStr = new String(arr);
						if (wordDict.contains(newStr)) {
							queue.add(newStr);
							wordDict.remove(newStr); // After visit remove the word
						}
						arr[i] = temp;
					}
				}
			}
			ladderLen++;
		}
		return 0;
	}

	/*
	 * Word Ladder II: Find ladder length and display the transformation
	 * Find all shortest transformation sequence(s) from beginWord to endWord, such that:
	 *    Input: beginWord = "hit", endWord = "cog"; wordList = ["hot","dot","dog","lot","log","cog"]
	 *    Output: [["hit","hot","dot","dog","cog"], ["hit","hot","lot","log","cog"]]
	 */
	public List<List<String>> wordLadderII(String start, String end, List<String> dict) {
		List<List<String>> result = new ArrayList<List<String>>();

		LinkedList<WordNode> queue = new LinkedList<WordNode>();
		queue.add(new WordNode(start, 1, null));

		// dict.add(end);

		int minStep = 0;

		HashSet<String> visited = new HashSet<String>();
		HashSet<String> unvisited = new HashSet<String>();
		unvisited.addAll(dict);

		int preNumSteps = 0;

		while (!queue.isEmpty()) {
			WordNode top = queue.remove();
			String word = top.word;
			int currNumSteps = top.count;

			if (word.equals(end)) {
				if (minStep == 0) minStep = top.count;

				if (top.count == minStep && minStep != 0) {
					// nothing
					ArrayList<String> t = new ArrayList<String>();
					t.add(top.word);
					while (top.prev != null) {
						t.add(0, top.prev.word);
						top = top.prev;
					}
					result.add(t);
					continue;
				}

			}

			// Why this???
			if (preNumSteps < currNumSteps) unvisited.removeAll(visited);

			preNumSteps = currNumSteps;

			char[] arr = word.toCharArray();
			for (int i = 0; i < arr.length; i++) {
				for (char c = 'a'; c <= 'z'; c++) {
					char temp = arr[i];
					if (arr[i] != c) {
						arr[i] = c;
					}

					String newWord = new String(arr);
					if (unvisited.contains(newWord)) {
						queue.add(new WordNode(newWord, top.count + 1, top));
						visited.add(newWord);
					}

					arr[i] = temp;
				}
			}
		}

		return result;
	}

	/* Unique Word Abbreviation:
	 * Assume you have a dictionary and given a word, find whether its abbreviation is unique in the dictionary. A word's abbreviation 
	 * is unique if no other word from the dictionary has the same abbreviation.
	 * Example:Given dictionary = [ "deer", "door", "cake", "card" ]
	 * isUnique("dear") -> false; isUnique("cart") -> true
	 */
	public boolean isUnique(Map<String, String> map, String word) {
		String abb = getAbbrevation(word);
		return (!map.containsKey(abb) || map.get(abb).equals(word));
	}

	public String getAbbrevation(String word) {
		int n = word.length();
		if (n <= 2) return word;

		return String.valueOf(word.charAt(0) + Integer.toString(n - 2) + word.charAt(n - 1));
	}

	public Map<String, String> buildDictionary(String[] strings) {
		Map<String, String> map = new HashMap<>();

		for (String str : strings) {
			String abb = getAbbrevation(str);
			// If there is any duplicate str, put empty or put str
			if (map.containsKey(abb)) map.put(abb, "");
			else map.put(abb, str);
		}

		return map;
	}

	/* Minimum Unique Word Abbreviation:
	 * A string such as "word" contains the following abbreviations: ["word", "1ord", "w1rd", "wo1d", "wor1", "2rd",
	 * "w2d", "wo2", "1o1d", "1or1", "w1r1", "1o2", "2r1", "3d", "w3", "4"] 
	 * Given a target string and a set of strings in a dictionary, find an abbreviation of this target string with the 
	 * smallest possible length such that it does not conflict with abbreviations of the strings in the dictionary. 
	 * Each number or letter in the abbreviation is considered length = 1. 
	 * For example, the abbreviation "a32bc" has length = 4. Examples: "apple", ["blade"] -> "a4" (because "5" or
	 * "4e" conflicts with "blade") "apple", ["plain", "amber", "blade"] -> "1p3" (other valid answers include "ap3", 
	 * "a3e", "2p2", "3le", "3l1").
	 */
	public String minAbbreviation(String target, String[] dictionary) {
		Set<String> visited = new HashSet<>();
		PriorityQueue<Abbr> q = new PriorityQueue<>((a, b) -> a.len - b.len);
		int len = target.length();
		String first = "";

		for (int i = 0; i < len; i++)
			first += "*";

		q.offer(new Abbr(first, 1));
		while (!q.isEmpty()) {
			Abbr ab = q.poll();
			String abbr = ab.abbr;
			boolean conflict = false;
			for (String word : dictionary) {
				if (word.length() == len && isConflict(abbr, word)) {
					conflict = true;
					break;
				}
			}
			if (conflict) generateAbbr(target, abbr, visited, q);
			else return NumAbbr(abbr);
		}

		return null;
	}

	boolean isConflict(String abbr, String str) {
		for (int i = 0; i < abbr.length(); i++)
			if (abbr.charAt(i) != '*' && str.charAt(i) != abbr.charAt(i)) return false;
		return true;
	}

	void generateAbbr(String str, String abbr, Set<String> visited, PriorityQueue<Abbr> q) {
		char[] temp = abbr.toCharArray();
		for (int i = 0; i < temp.length; i++) {
			if (temp[i] == '*') {
				temp[i] = str.charAt(i);
				String next = new String(temp);
				if (!visited.contains(next)) {
					q.offer(new Abbr(next, abbrLength(next)));
					visited.add(next);
				}
				temp[i] = '*';
			}
		}
	}

	int abbrLength(String abbr) {
		int ret = 0, star = 0;
		for (char c : abbr.toCharArray()) {
			if (c >= 'a' && c <= 'z') {
				ret += 1 + star;
				star = 0;
			} else if (c == '*') {
				star = 1;
			}
		}
		return ret + star;
	}

	String NumAbbr(String abbr) {
		String ret = "";
		int count = 0;
		for (char c : abbr.toCharArray()) {
			if (c >= 'a' && c <= 'z') {
				if (count > 0) {
					ret += count;
					count = 0;
				}
				ret += c;
			} else {
				count++;
			}
		}
		if (count > 0) ret += count;
		return ret;
	}

	/* Maximum Product of Word Lengths:
	 * Given a string array words, find the maximum value of length(word[i]) * length(word[j]) where the two words do not 
	 * share common letters. You may assume that each word will contain only lower case letters. If no such two words exist, return 0.
	 * Example 1: Input: ["abcw","baz","foo","bar","xtfn","abcdef"]	Output: 16
	 * Explanation: The two words can be "abcw", "xtfn".
	 * Example 2: Input: ["a","ab","abc","d","cd","bcd","abcd"]	Output: 4
	 * Explanation: The two words can be "ab", "cd".
	 * Example 3: Input: ["a","aa","aaa","aaaa"] Output: 0
	 * Explanation: No such pair of words.
	 * The solution is calculated by doing a product of the length of each string to every other string. Anyhow the
	 * constraint given is that the two strings should not have any common character. This is taken care by creating a
	 * unique number for every string. Image a an 32 bit integer where 0 bit corresponds to 'a', 1st bit corresponds to
	 * 'b' and so on. Thus if two strings contain the same character when we do and "AND" the result will not be zero
	 * and we can ignore that case.
	 */
	public int maxProduct(String[] words) {
		int n = words.length;
		int[] checker = new int[n];
		int masks = 0; // Variable to enable the bit based on the char
		for (int i = 0; i < n; i++) {
			masks = 0;
			for (int j = 0; j < words[i].length(); j++) {
				masks |= 1 << (words[i].charAt(j) - 'a');
			}
			checker[i] = masks;
		}

		int max = 0;
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				if ((checker[i] & checker[j]) == 0) // 0 means -> There is no same char b/w the strings
					max = Math.max(max, (words[i].length() * words[j].length()));
			}
		}

		return max;
	}

	/*
	 * Top K Frequent Words:
	 * Given a non-empty list of words, return the k most frequent elements. Your answer should be sorted by frequency
	 * from highest to lowest. If two words have the same frequency, then the word with the lower alphabetical order comes first.
	 * Example 1: Input: ["i", "love", "leetcode", "i", "love", "coding"], k = 2; Output: ["i", "love"]
	 */
	public void topKFrequentWords(String[] words, int k) {
		kElementsPattern.topKFrequentWords(words, k);
	}

	/*
	* Text Justification/Word Wrap Problem - GeeksforGeeks:
	* Given a sequence of words, and a limit on the number of characters that can be put in one line (line width). Put line breaks in the
	* given sequence such that the lines are printed neatly. Assume that the length of each word is smaller than the line width. 
	* The word processors like MS Word do task of placing line breaks. The idea is to have balanced lines. In other words, not have few 
	* lines with lots of extra spaces and some lines with small amount of extra spaces.
	* Eg: For example, consider the following string and line width M = 15
	*        "Geeks for Geeks presents word wrap problem"
	*     Following is the optimized arrangement of words in 3 lines
	*     		Geeks for Geeks
	*     		presents word
	*     		wrap problem 
	*/
	/*
	 * Approach1: The greedy solution is to place as many words as possible in the first line. Then do the same thing for the second line 
	 * and so on until all words are placed. This solution gives optimal solution for many cases, but doesn’t give optimal solution in all cases. 
	 * Despite being sub-optimal in some cases, the greedy approach is used by many word processors like MS Word and OpenOffice.org Writer.
	 * 
	 * Approach2: Recursive & DP
	 * The problem is to minimize the following total cost.
	 *  	Cost of a line = (Number of extra spaces in the line)^3 or (Number of extra spaces in the line)^2
	 *  	Total Cost = Sum of costs for all lines
	 * Please note that the total cost function is not sum of extra spaces, but sum of cubes (or square is also used) of extra spaces. 
	 * The idea behind this cost function is to balance the spaces among lines. For example, consider the following two arrangement of same 
	 * set of words: 
	 * 	1. There are 3 lines. One line has 3 extra spaces and all other lines have 0 extra spaces. Total extra spaces = 3 + 0 + 0 = 3. 
	 * 	   Total cost = 3*3*3 + 0*0*0 + 0*0*0 = 27.
	 * 	2. There are 3 lines. Each of the 3 lines has one extra space. Total extra spaces = 1 + 1 + 1 = 3. Total cost = 1*1*1 + 1*1*1 + 1*1*1 = 3.
	 * Ref: https://www.geeksforgeeks.org/word-wrap-problem-dp-19/
	 */

	//Using Greedy Algorithm:
	public String justify1(String words[], int width) {
		StringBuilder result = new StringBuilder();

		int i = 0, j = 0, n = words.length;
		while (i < n) {
			int len = words[i].length();
			result.append(words[i]);
			j = i + 1;
			//Here j-i-1 denotes no of space between words
			while (j < n && len + words[j].length() + (j - i - 1) < width) {
				result.append(" ").append(words[j]);
				len += words[j].length();
				j++;
			}
			if (j != n) result.append(" ");
			result.append("\n");
			i = j;
		}

		return result.toString();
	}

	/*
	 * Text Justification - Leetcode: Modification of above problem
	 * Given an array of words and a width maxWidth, format the text such that each line has exactly maxWidth characters 
	 * and is fully (left and right) justified.
	 * You should pack your words in a greedy approach; that is, pack as many words as you can in each line. Pad extra 
	 * spaces ' ' when necessary so that each line has exactly maxWidth characters.
	 */
	//Time: O(n*maxWidth); Space: O(lines*maxWidth); where n - no of words; lines - no of lines in the result
	public List<String> fullJustify(String[] words, int maxWidth) {
		int l = 0, r = 0, n = words.length;
		List<String> result = new ArrayList<>();
		while (l < n) {
			//1.Find the right pointer of words based on maxWidth
			int len = words[l].length();
			r = l + 1;
			//Here r-l-1 denotes no of space between words
			while (r < n && (len + words[r].length() + (r - l - 1)) < maxWidth) {
				len += words[r].length();
				r++;
			}

			//2.Justify the words
			result.add(justify(words, l, r - 1, len, maxWidth));
			l = r;
		}
		return result;
	}

	private String justify(String[] words, int l, int r, int wordsLength, int maxWidth) {
		//If it is single word in a line
		if (r - l == 0) return padResult(words[l], maxWidth);

		boolean isLastLine = r == words.length - 1;
		int noOfWords = r - l;
		int remSpace = maxWidth - wordsLength; //Remaining Space

		String space = isLastLine ? " " : whitespace(remSpace / noOfWords);
		int extraSpace = isLastLine ? 0 : remSpace % noOfWords;

		StringBuilder result = new StringBuilder();
		for (int i = l; i <= r; i++) {
			result.append(words[i]).append(space).append(extraSpace-- > 0 ? " " : "");
		}

		return padResult(result.toString().trim(), maxWidth);
	}

	private String padResult(String word, int maxWidth) {
		return word + whitespace(maxWidth - word.length());
	}

	public String whitespace(int numSpaces) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < numSpaces; i++) {
			sb.append(" ");
		}
		return sb.toString();
	}

	//Using DP:
	public String justify2(String words[], int width) {
		int cost[][] = new int[words.length][words.length];

		//Populate cost array
		populateCost(words, cost, width);

		//minCost from i to len is found by trying j between i to len and checking which one has min value
		int minCost[] = new int[words.length];
		int result[] = new int[words.length];
		for (int i = words.length - 1; i >= 0; i--) {
			minCost[i] = cost[i][words.length - 1];
			result[i] = words.length;
			for (int j = words.length - 1; j > i; j--) {
				if (cost[i][j - 1] == Integer.MAX_VALUE) continue;

				if (minCost[i] > minCost[j] + cost[i][j - 1]) {
					minCost[i] = minCost[j] + cost[i][j - 1];
					result[i] = j;
				}
			}
		}
		System.out.println("Minimum cost is " + minCost[0]);
		return getResult(words, result);
	}

	private void populateCost(String[] words, int[][] cost, int maxWidth) {
		int n = words.length;
		//Approach 1: using two for loops to calculate cost of putting words from i to j in one line. If words don't fit in one line then we put
		//Integer.MAX_VALUE there.
		for (int i = 0; i < n; i++) {
			cost[i][i] = maxWidth - words[i].length();
			for (int j = i + 1; j < n; j++) {
				cost[i][j] = cost[i][j - 1] - words[j].length() - 1;
			}
		}

		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				cost[i][j] = cost[i][j] < 0 ? Integer.MAX_VALUE : (int) Math.pow(cost[i][j], 2);
			}
		}

		// or 

		//Approach2: Combining two for loops into one 
		for (int i = 0; i < n; i++) {
			int width = maxWidth;
			for (int j = i; j < n; j++) {
				int len = width - words[j].length();
				cost[i][j] = len < 0 ? Integer.MAX_VALUE : (int) Math.pow(len, 2);
				width--; // For single space between words
			}
		}
	}

	private String getResult(String[] words, int[] result) {
		int i = 0;
		int j;
		//finally put all words with new line added in string buffer and print it.
		StringBuilder builder = new StringBuilder();
		do {
			j = result[i];
			for (int k = i; k < j; k++) {
				builder.append(words[k] + " ");
			}
			builder.append("\n");
			i = j;
		} while (j < words.length);
		return builder.toString();
	}

	// Crossword Puzzle: DFS & Backtracking
	public char[][] solvePuzzle(char[][] grid, String words) {
		return search(grid, Arrays.stream(words.split(";")).collect(Collectors.toSet()), 0, 0, 0);
	}

	static final int SIZE = 10;
	static final int[] R_OFFSETS = { 0, 1 };
	static final int[] C_OFFSETS = { 1, 0 };

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

	//TODO: Move this problem to appropriate category
	/*
	 * String Deletion: Given a String and dictionary hashset, write a function to determine the minimum no
	 * of characters to delete to make a word.
	 * Eg: dict = {a, aa, aaa}; 
	 * 		query = abc, o/p: 2; query: aac, o/p: 1;
	 */
	// Solution: BFS; Time-O(n!) because substring Time: n * n-1 * n-2... 1)
	public int stringDeletion(String query, HashSet<String> dict) {
		Queue<String> queue = new LinkedList<>();
		Set<String> visited = new HashSet<>();
		queue.add(query);
		while (!queue.isEmpty()) {
			String top = queue.poll();
			if (dict.contains(top)) return query.length() - top.length();
			// Check for all the substring
			for (int i = 0; i < top.length(); i++) {
				String subStr = top.substring(0, i) + top.substring(i + 1, top.length());
				if (subStr.length() > 0 && !visited.contains(subStr)) {
					visited.add(subStr);
					queue.add(subStr);
				}
			}
		}
		return -1;
	}
}

class Abbr {
	String abbr;
	int len;

	Abbr(String abbr, int len) {
		this.abbr = abbr;
		this.len = len;
	}
}