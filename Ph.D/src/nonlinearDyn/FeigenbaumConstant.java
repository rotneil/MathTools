package nonlinearDyn;

public class FeigenbaumConstant {
	// instance variable
	private static final int LOGISTIC = 1;
	private static final int HENON = 2;
	
	private double epsilon = 0.00000001;
	private int type = LOGISTIC;
	
	// interactive constructor
	public FeigenbaumConstant () {
		System.out.println("This program solves for the Feigenbaum constant in Logistic map and " +
				"Henon map");
		
		java.util.Scanner scanner = new java.util.Scanner(System.in);
		
		// prompt for map type
		System.out.println("Enter the type of the map (LOGISTIC MAP = 1, HENON MAP = 2): ");
		try {
			type = Integer.parseInt(scanner.next());
		} catch (NumberFormatException e) {
			System.out.println(e);
		}
		
		// allow multiple period computation
		int nextInput = 0;
		while (nextInput != -1) {
			// prompt for search parameters
			try {
				// prompt for the period
				System.out.print("Enter Period: ");
				int period = Integer.parseInt(scanner.next());
				
				// restrict period
				if (period > 1000)
					throw new IllegalArgumentException ("The entered period is beyond range!");
				
				// prompt for the lower limit
				System.out.print("Enter a1: ");
				double a1 = Double.parseDouble(scanner.next());
				
				
				// prompt for the upper limit
				System.out.print("Enter a2: ");
				double a2 = Double.parseDouble(scanner.next());
				System.out.println();
				
				double bPoint = getBifurcationPoint (period, a1, a2);
				System.out.println("Bifurcation Point: " + bPoint);
				
			} catch (IllegalArgumentException e) {
				System.out.println(e.toString());
			}
				
			System.out.println("Compute next bifurcation point. To stop iteration, enter -1");
			nextInput = Integer.parseInt(scanner.next());
			
			System.out.println("Enter the type of the map (LOGISTIC MAP = 1, HENON MAP = 2): ");
			type = Integer.parseInt(scanner.next());
			System.out.println();
		}
	}
	
	/**
	 * This method is a recursive method and it is used to find
	 * the bifurcation point for a certain period p withing the limit [a1, a2]
	 * @param p the period
	 * @param a1 the lower limit
	 * @param a2 the upper limit
	 * @return the bifurcation point
	 * @throws if the limit does not enclose the period given
	 */
	private double getBifurcationPoint (int period, double a1, double a2) 
			throws IllegalArgumentException
	{
		// get the period at each of the limits
		int p1 = getPeriod (a1);
		int p2 = getPeriod (a2);
		
		// ensure proper limit values
		if (period < p1)
			throw new IllegalArgumentException ("The given period is lower than the period " +
					p1 + " of the lower limit");
		if (period > p2)
			throw new IllegalArgumentException ("The given period is greater than the period " +
					p2 + " of the upper limit");
		
		// take the mid-value (a1, a2) as the running a-value
		double a = (a1 + a2) / 2.0;
		
		// base call
		if (Math.abs(a1 - a2) <= epsilon || p1 == p2)
			return a;
		
		// test for the period at this a-value
		int periodCount = getPeriod (a);
		
		// shift lower limit
		if (periodCount < period)
			return getBifurcationPoint (period, a, a2);
		
		// shift upper limit
		else
			return getBifurcationPoint (period, a1, a);
	}
	
	// method to test the period at a particular value of x
	private int getPeriod (double a)
	{
		// define the period counter
		int periodCount = 1;
		
		// distinguish between the two maps
		if (type == LOGISTIC) {
			// starting value of x
			double x = Math.random();
			
			// remove transient trajectory
			for (int i = 0; i < 20000000; ++i)
				x = f (a, x);
			
			// track the orbit to know the period
			double xx = x;
			x = f (a, x);
			while (Math.abs(xx - x) > 0.00000001) {
				x = f (a, x);
				++periodCount;
			}
			
		}
		else {
			// starting value of x
			double [] v = new double [] {Math.random(), Math.random()};
			
			// remove transient trajectory
			for (int i = 0; i < 25000000; ++i)
				v = f (a, v);
			
			// track the orbit to know the period
			double [] vv = new double [] {v[0], v[1]};
			v = f (a, v);
			while (Math.abs(vv[0] - v[0]) > 0.0000001) {
				v = f (a, v);
				++periodCount;
			}
			
		}
		
		return periodCount;
	}
	
	// method to define the logistic map function
	private double f (double a, double x)
	{
		// logistic map
		return (a * x * (1 - x));
	}
	
	// method to define Henon map
	private double [] f (double a, double [] v)
	{
		double [] vv = new double [v.length];
		vv [0] = a - v[0] * v[0] + 0.3 * v[1];
		vv [1] = v[0];
		
		return vv;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new FeigenbaumConstant ();
	}
}
