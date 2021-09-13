/**
 * This class contains the necessary implementation of a List
 * for use in the symbolic differentiation and other applications
 * like Postfix evaluation of expressions.
 * 
 * @author Nehemiah Oluwafemi
 * rotneil@yahoo.com
 */
package math.tools;


/**
 * A class to represent a node in the Lisk Data structure
 * @param <T> The data embedded in the list structure. Later, this
 * will be modified to handle operators and functions like tanh (x)
 */
class ListNode <T>
{
	// package access members
	T data;					// the data at the node
	ListNode <T> nextNode;	// reference to the next node in list
	
	/**
	 * The constructor creates a list node that has a generic object.
	 * This is used for link with just a data and no other listnode
	 * @param object The generic data at the node
	 */
	ListNode (T object)
	{
		this (object, null);
	}
	
	/**
	 * Constructor to create a ListNode that refers to the specified
	 * object and the next list
	 * @param object The data at the node
	 * @param node Reference to the next node
	 */
	ListNode (T object, ListNode <T> node)
	{
		data = object;
		nextNode = node;
	}
	
	/**
	 * Method to return the generic data at the node
	 * @return
	 */
	T getData ()
	{
		return data;
	}
	
	/**
	 * Method to return the reference to the next node
	 * @return
	 */
	ListNode <T> getNext ()
	{
		return nextNode;
	}
}	// end class ListNode<T>


/**
 * This is representation of data structures which has the ability
 * to increase at implementation.
 * @author Nehemiah Oluwafemi
 *
 * @param <T>
 */
public class List <T>
{
	// instance variables of the list
	private ListNode <T> firstNode;
	private ListNode <T> lastNode;
	private String name; 	// String for various impl of list
	
	/**
	 * No-argument constructor to give a default List 
	 */
	public List () {
		this ("List");
	}
	
	/**
	 * Argument constructor with just the name of the implementation
	 * @param listName
	 */
	public List(String listName) {
		this.name = listName;
		firstNode = lastNode = null;
	}
	
	/**
	 * This method inserts into the list from the front of the 
	 * list
	 * @param item
	 */
	public void insertAtFront (T item)
	{
		if (isEmpty ())
			firstNode = lastNode = new ListNode <T> (item);
		else
			firstNode = new ListNode<T> (item, firstNode);
	}
	
	/**
	 * This method inserts the argument from the back of the list
	 * @param item
	 */
	public void insertAtBack (T item)
	{
		if (isEmpty())
			firstNode = lastNode = new ListNode <T> (item);
		else
			lastNode = lastNode.nextNode = new ListNode<T> (item);
	}
	
	/**
	 * This method inserts an element into the list at the specified index.
	 * Note that indexing is zero-based.
	 * @param index
	 * @param item
	 * @throws IndexOutOfBoundsException If the size of List is less than index
	 * @throws EmptyListException
	 */
	public void insertAt (int index, T item)
		throws IndexOutOfBoundsException, EmptyListException
	{
		if (size () < index)
			throw new IndexOutOfBoundsException (String.format(
					"List size %d is less than index %d", size(), index));
		if (isEmpty())
			throw new EmptyListException (name + " is empty");
		
		// insert item at the index
		if (index == 0)					// At the front
			insertAtFront (item);
		else if (index == size())	// At the back
			insertAtBack (item);
		else {							// within the list
			int place = 0;
			ListNode <T> datum = new ListNode<T> (item);
			ListNode <T> current = firstNode;
			
			// loop until an element before the insertion point
			while (++place < index)
				current = current.nextNode;
			
			// break list at the point and insert element
			ListNode <T> nextList = current.nextNode;
			current.nextNode = datum;
			datum.nextNode = nextList;
		}
	}
	
	/**
	 * This method removes the first item from the list and returns
	 * the first element
	 * @return A reference to the removed first element
	 * @throws EmptyListException If the list is empty
	 */
	public T removeFromFront () throws EmptyListException
	{
		if (isEmpty())
			throw new EmptyListException (
					String.format("List %s is empty", name));
		
		// declare the first element
		T removedItem = firstNode.data;
		
		// update firstNode references
		if (firstNode == lastNode)
			firstNode = lastNode = null;
		else
			firstNode = firstNode.nextNode;
		
		return removedItem;
	}
	
	/**
	 * This method deletes the last element from the list and 
	 * returns the element
	 * @return A reference to the removed element
	 * @throws EmptyListException
	 */
	public T removeFromBack () throws EmptyListException
	{
		if (isEmpty ())
			throw new EmptyListException (
					String.format("List %s is empty", name));
		
		// get a reference to the data of the last node
		T removedItem = lastNode.data;
		
		// update the last node references
		if (firstNode == lastNode)
			firstNode = lastNode = null;
		else {
			// locate the new last node
			ListNode <T> current = firstNode;
			
			while (current.nextNode != lastNode)
				current = current.nextNode;
			
			lastNode = current;
			current.nextNode = null;
		}
		
		return removedItem;
	}
	
	/**
	 * This method removes an element at the specified index
	 * @param index
	 * @return
	 * @throws IndexOutOfBoundsException
	 * @throws EmptyListException
	 */
	public T removeElemenAt (int index) 
			throws IndexOutOfBoundsException, EmptyListException
	{
		if (size () < index)
			throw new IndexOutOfBoundsException (String.format(
					"List size %d is less than index %d", size(), index));
		if (isEmpty())
			throw new EmptyListException (name + " is empty");
		
		// insert item at the index
		if (index == 0)					// At the front
			return removeFromFront ();
		
		// At the back
		if (index == size())	
			return removeFromBack ();
		
		// within the list
		int place = 0;
		ListNode <T> current = firstNode;
		
		// loop until an element before the insertion point
		while (++place < index)
			current = current.nextNode;
		
		// break list at the point and insert element
		T element = current.nextNode.data;
		current = current.nextNode;
		return element;
	}
	
	/**
	 * This method returns a reference to the element at a 
	 * particular index without removing it from the list
	 * @param index
	 * @return
	 * @throws EmptyListException
	 */
	public T getElementAt (int index) throws EmptyListException,
		IndexOutOfBoundsException
	{
		// check for empty List
		if (isEmpty())
			throw new EmptyListException (
					String.format("Empty %s", name));
		if (index >= size())
			throw new IndexOutOfBoundsException (
					String.format("Elements in %s is less than" +
							" %d", name, index));
		ListNode <T> current = firstNode;
		
		int count = 0;
		while (++count < index)
			current = firstNode.nextNode;
		
		return current.data;
	}
	
	/**
	 * This method inserts into a List from the back. For instance, if the list already 
	 * has four elements, then the newly added element is the fifth.
	 * @param element
	 */
	public void add (T element)
	{
		insertAtBack (element);
	}
	
	/**
	 * Method that determines if a List is empty or not
	 * @return
	 */
	public boolean isEmpty () {
		return firstNode == null;
	}
	
	/**
	 * This is method returns the number of elements in the list
	 * @return
	 */
	public int size ()
	{
		if (isEmpty())
			return 0;
		
		ListNode <T> current = firstNode;
		int count = 0;
		while (current != null) {
			current = current.nextNode;
			++count;
		}
		return count;
	}
	
	@Override
	public String toString ()
	{
		if (isEmpty())
			return String.format("Empty %s\n", name);
		
		String output = String.format("The %s is: ", name);
		ListNode <T> current = firstNode;
		
		while (current != null) {
			output += current.data;
			current = current.nextNode;
		}
		
		return output;
	}
}
