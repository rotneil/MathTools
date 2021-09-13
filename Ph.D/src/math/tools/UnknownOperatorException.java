package math.tools;

public class UnknownOperatorException extends RuntimeException
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * No-argument constructor
	 */
	public UnknownOperatorException() {
		this ("Invalid operator");
	}
	
	public UnknownOperatorException (String name)
	{
		super (name);
	}
}
