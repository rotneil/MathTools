package math.tools;

public class EmptyListException extends RuntimeException
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * No-argument constructor
	 */
	public EmptyListException() {
		this ("List");
	}
	
	public EmptyListException (String name)
	{
		super (name);
	}
}
