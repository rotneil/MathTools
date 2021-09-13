/**
 * This is the implementation of Generic Data Structure Queue
 */

package math.tools;

public class Queue <T>
{
	// instance variable
	private List <T> queueList;
	
	/**
	 * No-argument constructor
	 */
	public Queue() {
		queueList = new List<T> ("Queue");
	}
	
	/**
	 * This method adds element to the queue from the back
	 * @param item
	 */
	public void enqueue (T item)
	{
		queueList.insertAtBack(item);
	}
	
	/**
	 * This method removes element from the queue from the front
	 * @return
	 * @throws EmptyListException
	 */
	public T dequeue () throws EmptyListException
	{
		return queueList.removeFromFront();
	}
	
	/**
	 * Method to check the content of a lis
	 * @return
	 */
	public boolean isEmpty ()
	{
		return queueList.isEmpty();
	}
	
	@Override
	public String toString () {
		return queueList.toString();
	}
}
