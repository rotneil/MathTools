package henonMap;

public class HenonLyapunov 
{
	// instance variables
	
	// constructor
	/**
	 * This instantiates the computation of HenonLyapunov exponents
	 */
	public HenonLyapunov (double a, double b)
	{
		// local variables
		int n = 200000;
		double [] h = new double [2];
		double [] v = new double [] {0.5, 0.5};		// the initial point
		
		// perform a hundred iteration
		for (int i = 0; i < 100; ++i)
			v = f(a, b, v);
		
		double [] [] w = getIdentity (2);
		double [] [] z;
		
		// perform the 
		for (int i = 0; i < n; ++i) {
			// get the new iterate
			v = f(a, b, v);
			
			// get the ellipsoid
			z = getProduct (Df(a, b, v), w);
			
			// orthogonalize z
			w = getOrthogonalVector (z);
			
			// get the expansion and convert to lyapunov exponents
			double [] len = getModulus (w);
			for (int j = 0; j < len.length; ++j)
				h [j] += Math.log(len[j]);
			
			// normalize the vector
			normalize (w);
		}
		
		for (int j = 0; j < h.length; ++j)
			System.out.println("h" + (j + 1) + " = " + (h[j] / n));
	}
	
	
	
	/**
	 * This method return the orthogonal matrix from z
	 * @param z
	 * @return
	 */
	private double [] [] getOrthogonalVector (double [] [] z)
	{
		// local variable
		double [] [] y = new double [z.length] [z[0].length];
		
		// put the values of z in y
		for (int i = 0; i < z.length; ++i)
			for (int j = 0; j < z[i].length; ++j)
				y[i][j] = z[i][j];
		
		// compute the orthogonal basis for the columns of z
		for (int j = 1; j < z[0].length; ++j) {
			for (int k = 0; k < j; ++k) {
				// compute the scalar product of the basis
				double dotSum = 0.0;
				double modulus = 0.0;
				
				// perform the column dot product and basis modulus
				for (int l = 0; l < z[j].length; ++l) {
					dotSum += z[l][j] * y[l][k];
					modulus += y[l][k] * y[l][k];
				}
				
				// orthogonalize
				for (int l = 0; l < z[j].length; ++l)
					y[l][j] -= dotSum * y[l][k] / modulus;
			}
		}
		
		return y;
	}
	
	
	/**
	 * This method normalizes the bases of the matrix y
	 * @param y
	 * @return
	 */
	private void normalize (double [] [] y)
	{
		// local variables
		double [] h = getModulus (y);
		
		// now divide each column by the modulus
		for (int j = 0; j < y[0].length; ++j)
			for (int i = 0; i < y.length; ++i)
				y [i][j] = y [i] [j] / h[j];
	}
	
	/**
	 * This method returns the length of the bases of the matrix y
	 * @param y
	 * @return
	 */
	private double [] getModulus (double [] [] y)
	{
		// local variable
		double [] h = new double [y.length];
		
		for (int j = 0; j < h.length; ++j) {
			for (int i = 0; i < y.length; ++i)
				h[j] = h[j] + y[i][j] * y[i][j];
			
			// find the square root
			h[j] = Math.sqrt(h[j]);
		}
		
		return h;
	}
	
	
	/**
	 * This method computes the product of two matrices a and b
	 * @param a
	 * @param b
	 * @return
	 */
	private double [] [] getProduct (double [] [] a, double [] [] b)
	{
		// local variable
		double [] [] c = new double [a.length] [b.length];
		
		for (int i = 0; i < a.length; ++i)
			for (int j = 0; j < b.length; ++j)
				for (int k = 0; k < a[i].length; ++k)
					c [i] [j] = c[i] [j] + a [i] [k] * b [k] [j];
		
		return c;
	}
	
	/**
	 * The computation of the Jacobi matrix
	 * @param a The a-value from Henon Map
	 * @param b The b-value of Henon Map
	 * @param v The present space phase
	 * @return The jacobi matrix
	 */
	private double [] [] Df(double a, double b, double [] v)
	{
		// local variables
		double [] [] D = new double [v.length] [v.length];
		D[0][0] = -2.0 * v[0];		D[0][1] = b;
		D[1][0] = 1.0;				D[1][1] = 0.0;
		
		return D;
	}
	
	/**
	 * The computation of the orbit
	 * @param a The a-value from Henon Map
	 * @param b The b-value of Henon Map
	 * @param v The previous orbit
	 * @return	The new orbit
	 */
	private double [] f (double a, double b, double [] v)
	{
		// declare local variable
		double [] vv = new double [v.length];
		vv [0] = a - v[0] * v[0] + b * v[1];
		vv [1] = v[0];
		
		return vv;
	}
	
	
	/**
	 * This method returns the identity matrix for n X n matrix
	 * @param n
	 * @return
	 */
	private double [] [] getIdentity (int n)
	{
		// local variable
		double [] [] in = new double [n] [n];
		
		for (int i = 0; i < n; ++i)
			in [i] [i] = 1.0;
		
		return in;
	}
	
	
	// launch the application
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		// prompt for the values of a and b
		java.util.Scanner scanner = new java.util.Scanner (System.in);
		System.out.println("Enter the value of a:");
		double a = scanner.nextDouble();
		
		System.out.println ("Enter the value of b:");
		double b = scanner.nextDouble();
		scanner.close();
		
		new HenonLyapunov(a, b);
	}

}
