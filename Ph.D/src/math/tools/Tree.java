/**
 * TreeNode and Tree class declarations for a binary search tree.
 */
package math.tools;

/**
 * TreeNode definition that handles the activities at each node
 * @author Nehemiah Oluwafemi
 *
 * @param <T>
 */
class TreeNode <T extends Comparable<T>>
{
	// package access members
	TreeNode <T> leftNode;
	T data;
	TreeNode <T> rightNode;
	
	/**
	 * Constructor initializes data and makes this a leaf node
	 * @param nodeData
	 */
	public TreeNode (T nodeData)
	{
		data = nodeData;
		leftNode = rightNode = null;
	}
	
	/**
	 * This method inserts a new value into a TreeNode by locating
	 * insertion point and inserting a new node;
	 * 
	 * Note: It ignores duplicate values
	 * @param data
	 */
	public void insert (T data)
	{
		// insert in left subtree
		if (data.compareTo(this.data) < 0) {
			// insert new TreeNode
			if (leftNode == null)
				leftNode = new TreeNode<T> (data);
			else	// continue traversing left subtree recursively
				leftNode.insert(data);
		}
		
		// insert in right subtree
		else if (data.compareTo(this.data) > 0) {
			// insert new TreeNode
			if (rightNode == null)
				rightNode = new TreeNode <T> (data);
			else   // traverse right subtree recursively
				rightNode.insert(data);
		}
	}
}

/**
 * Class definition for Tree
 * @author Nehemiah Oluwafemi
 *
 */
public class Tree <T extends Comparable<T>>
{
	// instance variable
	private TreeNode <T> root;
	
	/**
	 * Constructor initialize an empty Tree of data
	 */
	public Tree() {
		root = null;
	}

	/**
	 * This method inserts a new node in the binary tree
	 * @param data
	 */
	public void insertNode (T data)
	{
		if (root == null)
			root = new TreeNode <T> (data);	// create root node
		else
			root.insert(data); // call TreeNode method insert
	}
	
	/**
	 * This method deletes a node from the tree using the data provided.
	 * The method does it by seaching for the nodes that contains the item specified and
	 * then replacing it with the largest element of the subtree which is smaller than 
	 * the data in the node.
	 * 
	 * @param item
	 * @throws TreeNodeNotFoundException if no node if found containing the item
	 */
	public void deleteNode (T item) throws TreeNodeNotFoundException
	{
		// search for the node that contains the item
		TreeNode<T> deleteNode = contains (item);
		if (deleteNode == null)
			throw new TreeNodeNotFoundException ("Item not found in the tree");
		
		// get the deleteNode's parent
		TreeNode<T> deleteParent = getParentNode (deleteNode.data);
		
		// check if the node is a leaf
		if (deleteNode.rightNode == null && deleteNode.leftNode == null)
			replaceNode (deleteParent, deleteNode.data, null);
		// check if the node only has leftNode
		else if (deleteNode.rightNode == null)
			replaceNode (deleteParent, deleteNode.data, deleteNode.leftNode);
		// check if the node has only rightnode
		else if (deleteNode.leftNode == null)
			replaceNode (deleteParent, deleteNode.data, deleteNode.rightNode);
		else {
			// get the replacement node and its parent
			TreeNode<T> replaceParentNode = deleteNode;
			TreeNode<T> replaceNode = deleteNode.leftNode;
			while (replaceNode.rightNode != null) {
				replaceParentNode = replaceNode;
				replaceNode = replaceParentNode.rightNode;
			}
			
			// check if the replacement node is a leaf
			if (replaceNode.leftNode == null)
				replaceNode (replaceParentNode, replaceNode.data, null);
			else
				replaceNode (replaceParentNode, replaceNode.data, replaceNode.leftNode);
			
			// put the replacement node at the deleted node
			replaceNode (deleteParent, deleteNode.data, replaceNode);
			
			// put the left and right nodes of the deleteNode into replace node
			replaceNode.leftNode = deleteNode.leftNode;
			replaceNode.rightNode = deleteNode.rightNode;
		}
		
		deleteNode = null;
	}
	
	// method to replace a node with its childNode from the parent
	private void replaceNode (TreeNode<T> parent, T item, TreeNode<T> childNode)
	{
		if (parent == null)
			root = childNode;
		else if (parent.data.compareTo(item) < 0)
			parent.rightNode = childNode;
		else
			parent.leftNode = childNode;
	}
	
	/**
	 * This method searches for the parent node of the specified node. It returns null if
	 * the node is the root node.
	 * @param node
	 * @return
	 */
	public TreeNode<T> getParentNode (T item)
	{
		// compare with root's data
		int c = root.data.compareTo(item);
		
		if (c == 0)	// the element is found in root, and root cannot have a parent
			return null;
		
		// compare data with the node's data
		if (c < 0)
			return parentHelper (root, root.rightNode, item);
		else
			return parentHelper (root, root.leftNode, item);
	}
	
	// recursive method to handle node's parent
	private TreeNode<T> parentHelper (TreeNode<T> parent, TreeNode<T> node, T item)
	{
		// compare node's left and right data
		int c = node.data.compareTo(item);
		if (c == 0)
			return parent;
		else if (c < 0)
			return parentHelper (node, node.rightNode, item);
		else
			return parentHelper (node, node.leftNode, item);
	}
	
	
	/**
	 * This method returns the tree level. A value of 0 is returned if the tree is comprised
	 * of the root alone with no children nodes.
	 * @return
	 * @throws EmptyTreeException
	 */
	public int getDepth () throws EmptyTreeException
	{
		if (root == null)
			throw new EmptyTreeException ();
		
		return Math.max(depthHelper (root.leftNode), depthHelper (root.rightNode)); 
	}
	
	// Tree depth helper that recursively visits nodes
	private int depthHelper (TreeNode <T> node)
	{
		if (node == null)
			return 0;
		
		return (1 + Math.max(depthHelper(node.leftNode), depthHelper(node.rightNode)));
	}
	
	/**
	 * This method searches a Binary Tree for an item and returns the reference
	 * to the node containing the item. It returns null if that item cannot be located
	 * @param item
	 * @return
	 */
	public TreeNode <T> contains (T item)
	{
		if (root == null)
			return null;
		
		// compare this item with the root data
		int c = root.data.compareTo(item);
		
		if (c == 0)
			return root;
		else if (c < 0)
			return searchHelper(root.rightNode, item);
		else
			return searchHelper(root.leftNode, item);
	}
	
	// this is the recursive method that searches the Tree for a particular item.
	private TreeNode <T> searchHelper (TreeNode <T>node, T item)
	{
		if (node == null)
			return null;
		
		// compare node's data with item
		int comp = node.data.compareTo(item);
		
		if (comp == 0)
			return node;
		else if (comp < 0)
			return searchHelper (node.rightNode, item);
		else
			return searchHelper (node.leftNode, item);
	}
	
	/**
	 * This method traverses the tree level by level starting from the root and
	 * outputting levels from left to right.
	 * @return
	 */
	public String levelOrderTraversal () {
		// for empty tree
		if (root == null)
			return "";
		
		// enqueue the root
		String output = "";
		Queue<TreeNode<T>> queue = new Queue<TreeNode<T>>();
		queue.enqueue(root);
		
		// print from level to level while queue is not empty
		while (!queue.isEmpty()) {
			// dequeue and print it's data
			TreeNode<T> node = queue.dequeue();
			output += String.format("%s ", node.data);
			
			// enqueue node's children
			if (node.leftNode != null)
				queue.enqueue(node.leftNode);
			if (node.rightNode != null)
				queue.enqueue(node.rightNode);
		}
		
		return output;
	}
	
	/**
	 * This method traverses the Tree in preorder
	 * @return
	 */
	public String preorderTraversal () {
		return preorderHelper (root);
	}
	
	// recursive method to perform preorder traversal
	private String preorderHelper (TreeNode <T> node)
	{
		if (node == null)
			return "";
		
		String output = String.format("%s ", node.data);
		output += preorderHelper (node.leftNode);
		output += preorderHelper (node.rightNode);
		
		return output;
	}
	
	/**
	 * Method to traverse tree in order
	 * @return
	 */
	public String inorderTraversal () {
		return inorderHelper (root);
	}
	
	// recursive method to perform inorder traversal
	private String inorderHelper (TreeNode<T> node) {
		if (node == null)
			return "";
		
		String output = inorderHelper (node.leftNode);
		output += String.format("%s ", node.data);
		output += inorderHelper (node.rightNode);
		
		return output;
	}
	
	/**
	 * Method to output postorder elements
	 * @return
	 */
	public String postorderTraversal ()
	{
		return postorderHelper (root);
	}
	
	// recursive method to perform postorder traversal
	private String postorderHelper (TreeNode <T> node)
	{
		if (node == null)
			return "";
		
		String output = postorderHelper (node.leftNode);
		output += postorderHelper (node.rightNode);
		output += String.format("%s ", node.data);
		
		return output;
	}
	
	/**
	 * This method prints a 2 dimensional view of the tree with the rightmost subtree
	 * at the top of the right column of the view, the root at the middle of the left 
	 * column and the leftmost subtree at the bottom of the right column. Each column
	 * is separated by 5 spaces.
	 * @return
	 */
	public String outputTree ()
	{
		// the return string
		StringBuilder stringBuilder = new StringBuilder ();
		outputTree (stringBuilder, root, 0);
		return stringBuilder.toString();
	}
	
	// recursive method to help in the printing
	private void outputTree (StringBuilder lineOutput, TreeNode<T> node, int totalSpaces)
	{
		// while node is not null
		while (node != null) {
			outputTree (lineOutput, node.rightNode, (totalSpaces + 5));
			
			// print the white spaces
			for (int i = 0; i < totalSpaces; ++i)
				lineOutput.append(" ");
			
			// output the value in the current node
			lineOutput.append(String.format("%s\r\n", node.data));
			totalSpaces += 5;
			
			// set reference to refer to node's left subtree
			node = node.leftNode;
		}
	}
}
