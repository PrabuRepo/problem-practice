package com.consolidated.problems.design;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.LinkedList;
import java.util.Queue;

public class MaxQueue {

}

/*
 * Time Complexity:
 * 		Insertion: O(n) worse case, O(1) amortized because we remove at most N elements for N insertions.
 * 		Deletion and Max lookup: O(1) 
 * Space Complexity: O(n) on the Max queue
 */
class MaxQueue1 {
	Queue<Integer> queue;
	Deque<Integer> deque;

	public MaxQueue1() {
		queue = new LinkedList<>();
		deque = new ArrayDeque<>();
	}

	//Insert: O(n) worse case, O(1) amortized because we remove at most N elements for N insertions.
	public void enqueue(int data) {

		while (!deque.isEmpty() && data > deque.getLast()) {
			deque.removeLast();
		}

		queue.add(data);
		deque.addLast(data);
	}

	public int dequeue() {
		if (queue.isEmpty()) return -1;
		if (queue.peek() == max()) deque.removeFirst();

		return queue.poll();
	}

	public int max() {
		return queue.isEmpty() ? -1 : deque.getFirst();
	}

}