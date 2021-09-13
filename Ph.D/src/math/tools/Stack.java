
/**
 * This is the representation of First-In-Last-Out (FILO) 
 * implementation.
 * @author Nehemiah Oluwafemi
 *
 */

package math.tools;

public class Stack <T>
{
	// instance variables
	private List <T> stackList;
	
	/**
	 * No-argument constructor
	 */
	public Stack() {
		stackList = new List<T> ("Stack");
	}
	
	/**
	 * This method adds data to the stack from the front
	 * @param data
	 */
	public void push (T data) {
		stackList.insertAtFront(data);
	}
	
	/**
	 * This method removes the last inserted item into the stack
	 * @return
	 * @throws EmptyListException
	 */
	public T pop () throws EmptyListException {
		return stackList.removeFromFront();
	}
	
	/**
	 * This method returns the top element of the stack without
	 * popping it
	 * @return
	 * @throws EmptyListException
	 */
	public T peek () throws EmptyListException {
		return stackList.getElementAt(0);
	}
	
	/**
	 * This method determines if stack is empty
	 * @return
	 */
	public boolean isEmpty () {
		return stackList.isEmpty();
	}
	
	/**
	 * This method returns the amount of elements that are stacked in the list
	 * @return
	 */
	public int size () {
		return stackList.size();
	}
	
	@Override
	public String toString () {
		return stackList.toString();
	}
}
