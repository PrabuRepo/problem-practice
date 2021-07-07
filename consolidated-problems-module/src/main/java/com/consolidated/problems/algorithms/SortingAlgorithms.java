package com.consolidated.problems.algorithms;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;
import java.util.TreeMap;

import com.common.model.Interval;
import com.common.utilities.Utils;

public class SortingAlgorithms {

	/************************* Type1: Basic Sorting Problems ************************************/
	/*
	 * Merge Sorted Array:Given two sorted integer arrays nums1 and nums2, merge nums2 into nums1 as one sorted array.
	 */
	// Simple approach
	public void merge(int[] nums1, int m, int[] nums2, int n) {
		int i = m - 1, j = n - 1, index = nums1.length - 1;

		while (i >= 0 && j >= 0) nums1[index--] = (nums1[i] > nums2[j]) ? nums1[i--] : nums2[j--];

		while (j >= 0) nums1[index--] = nums2[j--];

		// System.out.println("Result: " + Arrays.toString(nums1));
	}

	public void merge2(int[] nums1, int m, int[] nums2, int n) {
		if (nums1.length == 0 || nums2.length == 0) return;
		int size = nums1.length - 1, i1 = m - 1, i2 = n - 1;
		// Merge from max value or last index
		while (size >= 0 && i1 >= 0 && i2 >= 0) {
			if (nums1[i1] > nums2[i2]) {
				nums1[size--] = nums1[i1--];
			} else {
				nums1[size--] = nums2[i2--];
			}
		}

		while (i2 >= 0) nums1[size--] = nums2[i2--];

		// System.out.println(Arrays.toString(nums1));
	}

	/*Largest Number: 
	 * Given a list of non negative integers, arrange them such that they form the largest number.
	 * Example 1: Input: [10,2]; Output: "210"; 
	 * Example 2: Input: [3,30,34,5,9]; Output: "9534330"
	 */
	//Time:  O(nklogn), Space:O(n); where n is length of array and k is average length of String;
	//Then compare 2 strings will take O(k).
	public String largestNumber(int[] nums) {
		if (nums.length == 0) return "0";

		String[] arr = new String[nums.length];
		for (int i = 0; i < nums.length; i++)
			arr[i] = String.valueOf(nums[i]);

		Arrays.sort(arr, (a, b) -> (b + a).compareTo(a + b));

		//After sorting if first value is zero, then all elements will be zero.
		if (arr[0].charAt(0) == '0') return "0";

		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < arr.length; i++)
			sb.append(arr[i]);

		return sb.toString();
	}

	//TODO: Maximum Gap - Bucket Sort/Radix Sort
	//TODO: Sort Big File/External Sort - Study this

	/*********************** Type2: Rearrangement Problems **************************************/
	// Sort Colors/Sort an array of 0s, 1s and 2s
	// 1.Using count array - With additional space
	public int[] sort012Approach1(int[] a) {
		int[] count = new int[3];

		for (int i = 0; i < a.length; i++)
			count[a[i]]++;

		int index = 0;
		while (count[0]-- > 0) a[index++] = 0;

		while (count[1]-- > 0) a[index++] = 1;

		while (count[2]-- > 0) a[index++] = 2;

		return a;
	}

	// 2.Using 3-way Partition
	public int[] sort012Approach2(int a[]) {
		int arr_size = a.length;
		int lo = 0, hi = arr_size - 1, mid = 0;
		while (mid <= hi) {
			if (a[mid] == 0) {
				Utils.swap(a, lo, mid);
				lo++;
				mid++;
			} else if (a[mid] == 1) {
				mid++;
			} else if (a[mid] == 2) {
				Utils.swap(a, mid, hi);
				hi--;
			}
		}

		return a;
	}
	// Sort Transformed Array - Parabola Prob - 2 ptr approach

	// Wiggle Sort I/Convert array into Zig-Zag fashion:A[0]<=A[1]>=A[2]<=A[3]>=A[4]<=A[5]
	public void wiggleSort(int[] nums) {
		if (nums == null || nums.length <= 1) {
			return;
		}

		for (int i = 0; i < nums.length - 1; i++) {
			if (i % 2 == 0) {
				if (nums[i] > nums[i + 1]) {
					Utils.swap(nums, i, i + 1);
				}
			} else {
				if (nums[i] < nums[i + 1]) {
					Utils.swap(nums, i, i + 1);
				}
			}
		}
	}

	// Wiggle Sort II/Sort an array in wave form/Peaks and Valleys -> arr[0]>=arr[1]<=arr[2]>=arr[3]<=arr[4]....
	void sortInWave(int arr[], int n) {
		// Traverse all even elements
		for (int i = 0; i < n; i += 2) {
			// If current even element is smaller
			// than previous
			if (i > 0 && arr[i - 1] > arr[i]) {
				Utils.swap(arr, i - 1, i);
			}

			// If current even element is smaller
			// than next
			if (i < n - 1 && arr[i] < arr[i + 1]) {
				Utils.swap(arr, i, i + 1);
			}
		}
	}

	/*Move all negative numbers to beginning and positive to end:
	 *  1. Not maintaing order
	 *  2. Maintaing order 
	 */
	/*
	 * 1.Don't need to maintain the order
	 *  	An array contains both positive and negative numbers in random order. Rearrange the array elements so that all negative 
	 *  numbers appear before all positive numbers.
	 *  Eg: Input : -12, 11, -13, -5, 6, -7, 5, -3, -6; Output :-12 -7 -3 -13 -5  -6 11 6 5
	 *  
	 *  Solution: Use Quick Sort Partition
	 */
	public void moveNegNumFront1(int arr[]) {
		int i = 0, j = 0, n = arr.length;
		while (j < n) {
			if (arr[j] < 0) {
				Utils.swap(arr, i, j);
				i++;
			}
			j++;
		}
	}

	/*
	 * 2.Maintain the order of elements
	 *  	An array contains both positive and negative numbers in random order. Rearrange the array elements so that all negative 
	 *  numbers appear before all positive numbers.
	 *  Eg: Input : -12, 11, -13, -5, 6, -7, 5, -3, -6; Output :-12 -13 -5 -7 -3 -6 11 6 5
	 */
	/*
	 * Solution1: A simple solution is to use another array. We copy all elements of original array to new array. We then traverse 
	 * the new array and copy all negative and positive elements back in original array one by one. Time: O(n^2); Space O(n)
	 * Solution2: Using Insertion Sort; Time: O(n^2); Space O(1)
	 * Solution3: Using Merge Sort; Time: O(nlogn); Space O(1); Merge Sort takes O(n) space, but here use reverse alg to save space;
	 */
	// Solution2: Insertion Sort Alg:
	public void moveNegNumFront21(int[] arr) {
		int key, j, n = arr.length;

		for (int i = 1; i < n; i++) {
			key = arr[i];
			// if current element is positive do nothing
			if (key > 0) continue;

			/* if current element is negative, shift positive elements of arr[0..i-1], to one position to their right */
			j = i - 1;
			while (j >= 0 && arr[j] > 0) {
				arr[j + 1] = arr[j];
				j = j - 1;
			}
			arr[j + 1] = key;
		}
	}

	// Solution3:
	public void moveNegNumFront22(int[] arr) {
		divideGroups(arr, 0, arr.length - 1);
	}

	public void divideGroups(int[] arr, int low, int high) {
		if (low >= high) return;
		int mid = (low + high) / 2;
		divideGroups(arr, low, mid);
		divideGroups(arr, mid + 1, high);
		mergeGroup(arr, low, mid, high);
	}

	/* Reverse Merge: steps to convert [Ln Lp Rn Rp] to [Ln Rn Lp Rp] without using extra space:
	 * 1. Reverse Lp and Rn. We get [Lp] -> [Lp'] and [Rn] -> [Rn'];    -> [Ln Lp Rn Rp] -> [Ln Lp� Rn� Rp]
	 * 2. Reverse [Lp� Rn�]. We get [Rn Lp];  =>  [Ln Lp� Rn� Rp] -> [Ln Rn Lp Rp]
	 */
	public void mergeGroup(int[] arr, int low, int mid, int high) {
		int l = low;
		int r = mid + 1;
		while (l <= mid && arr[l] <= 0) l++;
		while (r <= high && arr[r] <= 0) r++;

		reverse(arr, l, mid);
		reverse(arr, mid + 1, r - 1);
		reverse(arr, l, r - 1);
	}

	public void reverse(int[] arr, int l, int h) {
		while (l < h) {
			Utils.swap(arr, l, h);
			l++;
			h--;
		}
	}

	// Rearrange positive and negative values in an array - Quick Sort Partition
	/* 1.Rearrange positive and negative values in an array - Quick Sort Partition:
	 *  1. Not maintaing order
	 *  2. Maintaing order 
	 */

	/* 1. Not maintaing order
	 * An array contains both positive and negative numbers in random order. Rearrange the array elements so that positive and
	 * negative numbers are placed alternatively. The order of the appearance of elements is not maintained with this approach
	 * Eg: Input: [-1, 2, -3, 4, 5, 6, -7, 8, 9], then the output should be [9, -7, 8, -3, 5, -1, 2, 4, 6]
	 */
	public int[] rearrangePosAndNegNumbers1(int[] a) {
		// 1. Move all the negative numbers to front side
		int i = 0, j = 0;
		while (j < a.length) {
			if (a[j] < 0) {
				Utils.swap(a, i, j);
				i++;
			}
			j++;
		}

		/* 2.Now all positive numbers are at end and negative numbers at the beginning of array. Increment the negative 
		 *   index by 2 and positive index by 1, i.e. swap every alternate negative number with next positive number.
		 */
		int posIndex = i, negIndex = 0;
		while (posIndex < a.length && negIndex < posIndex && a[negIndex] < 0) {
			Utils.swap(a, posIndex, negIndex);
			posIndex++;
			negIndex += 2;
		}
		return a;
	}

	// 2. To maintain the order:
	/* Solution:
	 * The idea is to process array from left to right. While processing, find the first out of place element in the
	 * remaining unprocessed array. An element is out of place if it is negative and at odd index, or it is positive and
	 * at even index. Once we find an out of place element, we find the first element after it with opposite sign. We
	 * right rotate the subarray between these two elements (including these two).
	 */
	void rightrotate(int arr[], int n, int outofplace, int cur) {
		int tmp = arr[cur];
		for (int i = cur; i > outofplace; i--)
			arr[i] = arr[i - 1];
		arr[outofplace] = tmp;
	}

	public void rearrangePosAndNegNumbers2(int arr[]) {
		int outofplace = -1, n = arr.length;

		for (int index = 0; index < n; index++) {
			if (outofplace >= 0) {
				if (((arr[index] >= 0) && (arr[outofplace] < 0)) || ((arr[index] < 0) && (arr[outofplace] >= 0))) {
					rightrotate(arr, n, outofplace, index);

					// the new out-of-place entry is now 2 steps ahead
					if (index - outofplace > 2) outofplace = outofplace + 2;
					else outofplace = -1;
				}
			}

			// if no entry has been flagged out-of-place
			if (outofplace == -1) {
				// check if current entry is out-of-place: odd/even
				// if (((arr[index] >= 0) && ((index & 0x01) == 0)) || ((arr[index] < 0) && (index & 0x01) == 1))
				if ((arr[index] >= 0 && index % 2 == 0) || (arr[index] < 0 && index % 2 == 1)) {
					outofplace = index;
				}
			}
		}
	}

	/* Alternative Sorting:
	 * Given an array of integers, print the array in such a way that the first element is first maximum and second
	 * element is first minimum and so on.
	 * Eg: Input: arr[] = {1, 6, 9, 4, 3, 7, 8, 2}; Output : 9 1 8 2 7 3 6 4
	 */
	public void alternateSort(int arr[], int n) {
		Arrays.sort(arr);

		// Printing the last element of array first and then first element and then second last element and then second
		// element and so on.
		int i = 0, j = n - 1;
		while (i < j) {
			System.out.print(arr[j--] + " ");
			System.out.print(arr[i++] + " ");
		}

		// If the total element in array is odd then print the last middle element.
		if (n % 2 != 0) System.out.print(arr[i]);
	}

	// Relative Sorting - Sorting based on another array

	/****************** Type3: Min no of swap required to sort array ***********************/
	//Approach1: Brute Force Approach: Time: O(n), Space: O(n)
	public int findUnsortedSubarray1(int[] nums) {
		int begin = nums.length - 1;
		int end = nums.length - 1;

		int[] sorted = nums.clone();
		Arrays.sort(sorted);

		for (int i = 0; i < nums.length; i++) {
			if (sorted[i] != nums[i]) {
				begin = i;
				break;
			}
		}
		for (int i = nums.length - 1; i > begin; i--) {
			if (sorted[i] != nums[i]) {
				end = i;
				break;
			}
		}

		return (end == begin) ? 0 : (end - begin + 1);
	}

	//Approach2: Time: O(n), Space: O(1)
	public int findUnsortedSubarray2(int[] nums) {
		int n = nums.length, start = 0, end = n - 1;

		//Find first dip from start
		while (start < n - 1) {
			if (nums[start] > nums[start + 1]) break;
			start++;
		}

		// There is no unsorted subarray
		if (start == n - 1) return 0;

		//Find last bump from end 
		while (end > 0) {
			if (nums[end] < nums[end - 1]) break;
			end--;
		}

		//Find the min and max value between start and end index
		int max = Integer.MIN_VALUE, min = Integer.MAX_VALUE;
		for (int i = start; i <= end; i++) {
			max = Math.max(max, nums[i]);
			min = Math.min(min, nums[i]);
		}

		//Expand start and end outward:
		//Find any other element greater than min in the left side
		while (start > 0 && nums[start - 1] > min) {
			start--;
		}

		//Find any other element less than max in the right side
		while (end < n - 1 && nums[end + 1] < max) {
			end++;
		}

		return end - start + 1;
	}

	//Approach3: Using Montonic Stack: Time: O(n), Space: O(1)
	//Best example to understand the monotonic stack
	public int findUnsortedSubarray3(int[] nums) {
		Stack<Integer> stack = new Stack<>();
		int l = nums.length, r = 0;
		//monotonous increasing stack - Find start index
		for (int i = 0; i < nums.length; i++) {
			while (!stack.isEmpty() && nums[stack.peek()] > nums[i]) l = Math.min(l, stack.pop());
			stack.push(i);
		}
		stack.clear();

		//monotonous decreasing stack - Find end index
		for (int i = nums.length - 1; i >= 0; i--) {
			while (!stack.isEmpty() && nums[stack.peek()] < nums[i]) r = Math.max(r, stack.pop());
			stack.push(i);
		}
		return r - l > 0 ? r - l + 1 : 0;
	}

	/********************************* Type4: Cyclic Sort/Marker Alg ************************/
	// 6.Find the Corrupt Pair
	/*
	 * Given an array containing n+1 numbers taken from the range 1 to n. One of the numbers got duplicated which also resulted in one number going
	 * missing. Find these numbers.
	 * Approach1: Cyclic Sort:
	 * Approach2: Marker Approach
	 */
	public int[] findCorruptPair(int[] nums) {
		// rearrange the array using cyclic sort.
		int i = 0, n = nums.length;
		while (i < n) {
			int val = nums[i] - 1;
			if (nums[val] != nums[i]) Utils.swap(nums, val, i);
			else i++;
		}

		for (i = 0; i < nums.length; i++)
			if (nums[i] != i + 1) return new int[] { nums[i], i + 1 };

		return new int[] { 0, 0 };
	}

	// 8.Find the First K Missing Positive Numbers (hard)

	// 9.Insert into a Cyclic Sorted List - Additional Prob - Check this

	/************************* Type: Revisit and Add Category ************************************/

	/*
	 * Given a sorted array in non-decreasing order, return an array of squares of each number, also 
	 * in non-decreasing order.
	 * For example: [-4,-2,-1,0,3,5] -> [0,1,4,9,16,25]
	 */
	/*
	 * Solution:
	 * Approach1: Sorting Approach: Time: O(nlogn), Space: O(1)
	 * 	Doing this problem in O(nlog(n)) time is pretty trivial - just square all the numbers and then
	 * 	sort them.
	 *
	 * Approach2: Linear Solution -> Time: O(n), Space: O(n)
	 *	There is a pattern in the input array. The largest squares will be at either end of the array.
	 *	The lowest -ve number and the highest +ve number will be the largest squares. So, if we look at
	 *	either ends of the array, we can go inwards and find smaller squares. This will give us squares 
	 *  in descending order - from largest to smallest.
	 *  Keep in mind that we will need to store the output somewhere. We will need to allocate a separate
	 *  array for that. Unfortunately, we cannot do this in-place (i.e, by rearranging the input array).
	 *  We allocate a new array and fill it from the back (since our squares are presented from largest to
	 *  smallest).
	 *	
	 */
	public int[] sortedSquares(int[] arr) {
		if (arr == null || arr.length <= 1) return arr;

		int n = arr.length, l = 0, h = n - 1, i = n - 1;
		int[] result = new int[n];

		while (l <= h) {
			if (Math.abs(arr[l]) >= Math.abs(arr[h])) {
				result[i--] = arr[l] * arr[l];
				l++;
			} else {
				result[i--] = arr[h] * arr[h];
				h--;
			}
		}

		return result;
	}

	/********************* Type1: Interval Patterns - Selection Problems **************************/
	/*
	 * Max length chain/Maximum Length of Pair Chain: 
	 * Time Complexity: O(nlogn)
	 */
	public int findLongestChain(int[][] pairs) {
		int count = 1, i = 0, j = 1;
		Arrays.sort(pairs, (a, b) -> a[0] - b[0]);
		while (i < pairs.length && j < pairs.length) {
			if (pairs[i][1] < pairs[j][0]) {
				count++;
				i = j;
				j++;
			} else {
				if (pairs[i][1] > pairs[j][1]) i = j;
				j++;
			}
		}

		return count;
	}

	/******************** Type1: Interval Patterns - Interval Manipulations ***********************/

	/* Data Stream as Disjoint Intervals:
	 * Given a data stream input of non-negative integers a1, a2, ..., an, ..., summarize the numbers seen so far as a list
	 * of disjoint intervals. For example, suppose the integers from the data stream are 1, 3, 7, 2, 6, ..., then the
	 * summary will be: 
	 * [1, 1] 
	 * [1, 1], [3, 3] 
	 * [1, 1], [3, 3], [7, 7] 
	 * [1, 3], [7, 7] 
	 * [1, 3], [6, 7]
	 */
	/** Initialize your data structure here. */
	LinkedList<Integer> list = new LinkedList<>();
	TreeMap<Integer, Interval> tree = new TreeMap<>();

	// Brute Force Approach
	public void addNum1(int val) {
		/*list.add(val);
		Collections.sort(list);*/
		int i = 0;
		for (i = 0; i < list.size(); i++) {
			if (val == list.get(i)) return;
			if (val < list.get(i)) {
				list.add(i, val);
				break;
			}
		}
		if (list.isEmpty() || list.size() == i) list.add(val);
	}

	public List<Interval> getIntervals1() {
		List<Interval> intervals = new ArrayList<>();
		for (int i = 0; i < list.size(); i++) {
			int s = i;
			while (i < list.size() - 1 && (list.get(i) + 1 == list.get(i + 1) || list.get(i) == list.get(i + 1))) i++;

			intervals.add(new Interval(list.get(s), list.get(i)));
		}
		return intervals;
	}

	public void addNum2(int val) {
		if (tree.containsKey(val)) return;
		Integer l = tree.lowerKey(val);
		Integer h = tree.higherKey(val);
		if (l != null && h != null && tree.get(l).end + 1 == val && h == val + 1) {
			tree.get(l).end = tree.get(h).end;
			tree.remove(h);
		} else if (l != null && tree.get(l).end + 1 >= val) {
			tree.get(l).end = Math.max(tree.get(l).end, val);
		} else if (h != null && h == val + 1) {
			tree.put(val, new Interval(val, tree.get(h).end));
			tree.remove(h);
		} else {
			tree.put(val, new Interval(val, val));
		}
	}

	public List<Interval> getIntervals2() {
		return new ArrayList<>(tree.values());
	}

	//TODO: Move below problems to appropriate category

	/*
	 * Triplet Sum: Given 3 arrays a,b,c of different sizes, find the number of distinct triplets(p,q,r) 
	 * where p is an element of a, written as p->a,q->b and r->c, satisfying the criteria: p<=q && q>=r.
	 * For example, a={3,5,7} b={3,6} and c={4,6,9}, we find four distinct triplets: (3,6,4), (3,6,6), (5,6,4), (5,6,6)
	 */
	// Approach1: Brute Force Approach
	static long triplets1(int[] a, int[] b, int[] c) {
		int count = 0, p, q, r;

		for (int i = 0; i < a.length; i++) {
			p = a[i];
			if (i > 0 && a[i - 1] == a[i]) continue;
			for (int j = 0; j < b.length; j++) {
				q = b[j];
				if (p > q || (j > 0 && b[j - 1] == b[j])) continue;
				for (int k = 0; k < c.length; k++) {
					r = c[k];
					if (k > 0 && a[k - 1] == a[k]) continue;

					if (q >= r) count++;
				}
			}
		}
		return count;
	}

	// Approach2: Sorting & compare with b[] array- Time Complexity-O(n^3)
	static long triplets(int[] a, int[] b, int[] c) {
		Arrays.sort(a);
		Arrays.sort(b);
		Arrays.sort(c);

		int p = 0, r = 0;
		long pCount = 0, rCount = 0, total = 0;
		for (int q = 0; q < b.length; q++) {
			while (p < a.length && a[p] <= b[q]) {
				if (p == 0 || a[p - 1] != a[p]) pCount++;
				p++;
			}
			while (r < c.length && c[r] <= b[q]) {
				if (r == 0 || c[r - 1] != c[r]) rCount++;
				r++;
			}
			if (q == 0 || b[q - 1] != b[q]) total += pCount * rCount;
		}

		return total;
	}

}