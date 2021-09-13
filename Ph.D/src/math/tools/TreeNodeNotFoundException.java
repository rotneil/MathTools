package math.tools;

public class TreeNodeNotFoundException extends RuntimeException
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * No-argument constructor
	 */
	public TreeNodeNotFoundException() {
		this ("TreeNode is not found.");
	}
	
	public TreeNodeNotFoundException (String name)
	{
		super (name);
	}
}
