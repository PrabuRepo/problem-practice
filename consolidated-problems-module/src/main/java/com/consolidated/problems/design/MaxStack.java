package com.consolidated.problems.design;

import java.util.Stack;

public class MaxStack {

	public static void main(String[] args) {
		MaxStack ob = new MaxStack();
		ob.testStackMax1();
		ob.testStackMax2();
	}

	public void testStackMax1() {
		MaxStack1 stack = new MaxStack1();
		System.out.println("Stack (using Auxiliary Stack): ");
		stack.push(4);
		stack.push(3);
		stack.push(5);
		stack.push(1);
		stack.push(2);
		System.out.println("Max: " + stack.getMax());
		System.out.println("Pop: " + stack.top());
		stack.pop();
		System.out.println("Max: " + stack.getMax());
		System.out.println("Pop: " + stack.top());
		stack.pop();
		System.out.println("Max: " + stack.getMax());
		System.out.println("Pop: " + stack.top());
		stack.pop();
		System.out.println("Max: " + stack.getMax());
		System.out.println("Pop: " + stack.top());
		stack.pop();
	}

	public void testStackMax2() {
		MaxStack21 stack = new MaxStack21();
		System.out.println("Stack (using Max variable): ");
		stack.push(4);
		stack.push(3);
		stack.push(5);
		stack.push(1);
		stack.push(2);
		System.out.println("Max: " + stack.max());
		System.out.println("Pop: " + stack.pop());
		System.out.println("Max: " + stack.max());
		System.out.println("Pop: " + stack.pop());
		System.out.println("Max: " + stack.max());
		System.out.println("Pop: " + stack.pop());
		System.out.println("Max: " + stack.max());
		System.out.println("Pop: " + stack.pop());
	}
}

// Using additional Min Variable
/*
 * Time Complexity: O(1) for all the operations
 * Space Complexity: O(n) 
 */
class MaxStack1 {
	public StckNode top;

	public MaxStack1() {

	}

	public void push(int x) {
		if (top == null) {
			top = new StckNode(x, x);
		} else {
			StckNode node = new StckNode(x, Math.max(x, top.max));
			node.next = top;
			top = node;
		}
	}

	public void pop() {
		if (top == null) return;
		top = top.next;
	}

	public int top() {
		if (top == null) return -1;
		return top.value;
	}

	public int getMax() {
		if (top == null) return -1;
		return top.max;
	}
}

class StckNode {
	public int value;
	public int max;
	public StckNode next;

	public StckNode(int value, int max) {
		this.value = value;
		this.max = max;
	}
}

// Two Stacks Approach: Design and Implement Special Stack Data Structure, using one auxiliary max stack
/*
 * Time Complexity: O(1) for all the operations
 * Space Complexity: O(n) for additional stack
 */
class MaxStack21 extends Stack<Integer> {
	private static final long serialVersionUID = 53463461L;

	// To hold the max values
	Stack<Integer> maxStack;

	public MaxStack21() {
		maxStack = new Stack<>();
	}

	public void push(int data) {
		if (this.isEmpty() || data >= max()) maxStack.push(data);
		// Push into main stack
		super.push(data);
	}

	public Integer pop() {
		if (super.isEmpty()) return -1;

		if (super.peek() == max()) maxStack.pop();
		return super.pop();
	}

	public int max() {
		return maxStack.isEmpty() ? Integer.MIN_VALUE : maxStack.peek();
	}
}

//Two Stacks Approach: Modification of above
class MaxStack22 {

	//To hold all the values
	Stack<Integer> stack;

	// To hold the max values
	Stack<Integer> maxStack;

	public MaxStack22() {
		stack = new Stack<>();
		maxStack = new Stack<>();
	}

	public void push(int data) {
		if (stack.isEmpty() || data >= max()) maxStack.push(data);
		// Push into main stack
		stack.push(data);
	}

	public Integer pop() {
		if (stack.isEmpty()) return -1;

		if (stack.peek() == max()) maxStack.pop();
		return stack.pop();
	}

	public int max() {
		return maxStack.isEmpty() ? Integer.MIN_VALUE : maxStack.peek();
	}
}