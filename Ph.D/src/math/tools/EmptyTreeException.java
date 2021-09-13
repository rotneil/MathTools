package math.tools;

public class EmptyTreeException extends RuntimeException
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * No-argument constructor
	 */
	public EmptyTreeException() {
		this ("Tree is empty");
	}
	
	public EmptyTreeException (String name)
	{
		super (name);
	}
}
