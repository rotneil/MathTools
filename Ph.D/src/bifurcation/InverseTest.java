package bifurcation;

public class InverseTest {

	public InverseTest()
	{
		// TODO Auto-generated constructor stub
		
		// obtain inputs for the matrix
		double [] [] a;
		
		java.util.Scanner scanner = new java.util.Scanner(System.in);
		try {
			System.out.println("Enter Matrix dimension: ");
			int dim = scanner.nextInt();
			
			while (dim != -1) {
				a = new double [dim] [dim];
				
				// accept inputs for the rows
				for (int i = 0; i < dim; ++i) {
					System.out.println("Enter the elements of the " + (i + 1) + " row");
					for (int j = 0; j < dim; ++j) {
						a[i][j] = scanner.nextDouble();
						System.out.println("Next element");
					}
				}
				System.out.print("\nMatrix A is: \n\t");
				displayMatrix (a);
				
				// evaluate inverse
				double [] [] inverse = inverse (a);
				System.out.print("\nThe inverse is: \n\t");
				
				displayMatrix(inverse);
				
				// prompt for the system dimension again
				System.out.println("Enter Matrix dimension: ");
				dim = scanner.nextInt();	
			}
		} catch (Exception e) {
			System.out.println (e.getMessage());
		}
	}

	private void displayMatrix (double [] [] a) {
		int dim = a.length;
		
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j < dim; ++j)
				System.out.print (a[i][j] + "    ");
				
			System.out.print("\n\t");
		}
		
		System.out.println();
	}
	/**
	 * This method calls the standard inverse method which uses the Gauss Elimination method
	 * for evaluation of inverse by setting the pivot and scale to default and the tolerance to
	 * 10(-12)
	 * @param a The matrix whose inverse is desired
	 * @return The inverse of the matrix a
	 * @throws MatrixException
	 */
	private double [] [] inverse (double [] [] a) throws MatrixException
	{
		// make a copy of a
		double [] [] aCopy = a.clone();
		
		// default values
		return inverse(aCopy, new int [a.length], new double [a.length], 0.000000000001); 
	}
	
	/**
	 * This method uses Gauss Elimination method to determine the inverse of a matrix [A] that 
	 * is a square matrix.
	 * @param A The matrix whose inverse is desired
	 * @param pivot The transformation sequence of the rows during pivot operation
	 * @param scale The scaling factor to enhance elimination process
	 * @param tol The tolerance, which should be as low as possible to avoid division by zero
	 * @return The inverse of matrix A i.e [A]'
	 * @throws MatrixException If any error occurs, especially if a value less than tol is obtained
	 */
	private double [] [] inverse (double [] [] a, int [] pivot, double [] scale, double tol) 
			throws MatrixException
	{
		// local variable that separates the argument from the return
		double [] [] x = new double [a.length][a.length];
		double [] b = new double [a.length];
		
		// decompose matrix a
		decompose(a, pivot, scale, tol);
		
		// Prepare the decomposed matrix for substitution process
		// Given [A]{x} = {B}, let the elements of matrix {x} represents the inverse of [A].
		// Then, those of {B} are of a unit matrix
		for (int i = 0; i < a.length; ++i) {
			for (int j = 0; j < b.length; ++j)
				b[j] = (i == j ? 1.0 : 0.0);
			
			// call substitution method
			substitute(a, pivot, b, x[i]);
		}
		
		// return the transposed value of x
		return transpose (x);
	}
	
	/**
	 * This method uses Gauss elimination method to seperate matrix A into the corresponding
	 * lower and upper triangular matrix. The elements of the formed upper triangular matrix
	 * represent the U(i,j) from i = 0 to n and j = i to n, while the lower triangular matrix
	 * are L(i, j) from i = 1 to n for j = 0 to n - 1. The diagonal element of L are all 1's.
	 * @param a The matrix to be decomposed into upper and lower matrices
	 * @param pivot The pivot transformation order. During decomposition, rows are supposed to
	 * be interchanged to enhance the elimination process. But instead, rows are retained but 
	 * the indices of the pivot are stored in this argument
	 * @param scale The scaling factor of the matrix. This is the maximum element in a row
	 * @param tol The tolerance, below which an element is assumed to be zero
	 * @return The decomposed matrix that has the upper and lower triangular matrix
	 * @throws MatrixException
	 */
	private double [] [] decompose (double [] [] a, int [] pivot, double [] scale, double tol)
			throws MatrixException
	{
		// confirm the validity of matrix a
		if (a.length != a[0].length)
			throw new MatrixException ("Matrix a is not a square matrix");
		if (pivot.length != a.length)
			throw new MatrixException ("The rows of matrix a must equal to the length of pivot");
		if (scale.length != a.length)
			throw new MatrixException ("The rows of matrix a must equal to the length of scale");
		
		// local variable
		double [] [] c = new double [a.length][a.length];
		
		// set the scale and the pivot to reflect the normal indices
		for (int i = 0; i < a.length; ++i) {
			pivot[i] = i;
			
			// set the first element of the row to scale
			scale [i] = a[i][0];
			
			// choose the maximum element in the row as the scale
			for (int j = 1; j < a.length; ++j)
				if (Math.abs(a[i][j]) > scale [i])
					scale [i] = a[i][j];
		}
		
		// now begin the elimination process
		for (int k = 0; k < a.length - 1; ++k) {
			// call method setPivot
			setPivot (a, pivot, scale, k);
			
			// check the system's tolerance before elimination
			double scaledValue = Math.abs(a[pivot[k]][k]) / scale[pivot[k]];
			if (scaledValue < tol)
				throw new MatrixException ("Scaled value " + scaledValue + 
						"is almost equal to zero: Less than the required tolerance " + tol);
			
			// operate on the lower rows
			for (int i = k + 1; i < a.length; ++i) {
				// evaluate the multiplication factor
				double factor = a [pivot[i]][k] / a [pivot[k]][k];
				
				// set the lower triangular matrix element
				a[pivot[i]][k] = factor;
				
				// eliminate pivot row from the current row
				for (int j = k + 1; j < a.length; ++j)
					a[pivot[i]][j] = a [pivot[i]][j] - factor * a [pivot[k]][j];
			}
		}
		
		return c;
	}
	
	/**
	 * This method sets the pivot for the Gauss elimination method.
	 * @param a This is the Matrix to be be decomposed by Gauss elimination method
	 * @param pivot This is the row transformation sequence that is to be set for this row operation
	 * @param scale The scaling value for the rows
	 * @param row The current Row operation that requires pivoting
	 */
	private void setPivot (double [] [] a, int [] pivot, double [] scale, int index)
	{
		// local variable
		double dummy;
		
		// set the initial value for pivotRow and set the max element
		int rowPivot = index;
		double max = Math.abs(a[pivot[index]][index]) / scale[pivot[index]];
		
		// iterate through the elements of this column for the biggest
		for (int i = index + 1; i < a.length; ++i) {
			// evaluate scaled value of the elements as dummy and compare it with other elements
			dummy = Math.abs(a[pivot[i]][index]) / scale[pivot[i]];
			if (dummy > max) {
				max = dummy;
				rowPivot = i; 
			}
		}
		
		// interchange the content of pivot at index and rowPivot
		int dummyIndex = pivot [rowPivot];
		pivot [rowPivot] = index;
		pivot [index] = dummyIndex;
	}	// end method setPivot
	
	/**
	 * This method implements the forward and backward substitution of a decomposed matrix A 
	 * in order to evaluate the elements of the unknown Matrix X where the 
	 * system of linear equation is given as
	 * 		[A]{X} = {B}
	 * The forward substitution implements
	 * 		[L]{D} = {B},
	 * while the backward then finds the unknown [X]
	 * 		[U]{x} = {D}
	 * The row transformation order is stored in matrix pivot.
	 * @param a The decomposed Upper and lower triangular matrix
	 * @param pivot The transformation sequence
	 * @param b The rigth-hand-side of the system of linear equation
	 * @param x The solution of the matrix
	 */
	private void substitute (double [] [] a, int [] pivot, double [] b, double [] x)
	{
		// local variable
		int n = a.length;
		double sum;
		
		// begin the forward substitution
		for (int i = 0; i < n; ++i) {
			sum = b[pivot[i]];
			
			// evaluate the summation and perform the subtraction
			for (int j = 0; j < i; ++j)
				sum = sum - a[pivot[i]][j] * b[pivot[j]];
			
			// set the value of d for i
			b[pivot[i]] = sum;
		}
		
		// begin the backward substitution
		x[n - 1] = b[pivot[n - 1]] / a [pivot[n - 1]][n - 1];
		
		// evaluate the values of matrix {x}
		for (int i = n - 2; i >= 0; --i) {
			sum = 0.0;
			
			// compute the summation
			for (int j = i + 1; j < n; ++j)
				sum = sum + a [pivot[i]][j] * x[j];
			
			x[i] = (b[pivot[i]] - sum) / a[pivot[i]][i];
		}
	}
	
	// method to perform multiplication of a n x n matrix by a row matrix
	private double [] [] product (double [] [] a, double [] b) throws MatrixException
	{
		// transpose matrix b and then call full product
		return product (a, transpose (b));
	}
	
	// method to perform full matrix multiplication
	private double [] [] product (double [] [] a, double [] [] b) throws MatrixException
	{
		// check validity of the arguments
		if (a == null || b == null)
			throw new MatrixException ("Matrix " + (a == null ? "A " : "B ") + "is empty");
		if (a[0].length != b.length)
			throw new MatrixException ("The columns of matrix A is not equal to rows of B");
		
		// local variables
		double [] [] c = new double [a.length] [b[0].length];
		
		for (int i = 0; i < a.length; ++i)
			for (int j = 0; j < b[i].length; ++j)
				for (int k = 0; k < a[i].length; ++k)
					c[i][j] = c[i][j] + a[i][k] * b [k][j];
		
		return c;
	}
	
	// method to transpose i x j matrix
	private double [] [] transpose (double [] [] a)
	{
		double [] [] t = new double [a[0].length] [a.length];
		
		for (int i = 0; i < a.length; ++i)
			for (int j = 0; j < a[i].length; ++j)
				t [i] [j] = a [j][i];
		
		return t;
	}
	
	// method to transpose a row matrix
	private double [] [] transpose (double [] a)
	{
		double [] [] t = new double [a.length][1];
		
		for (int i = 0; i < a.length; ++i)
			t [i][0] = a[i];
		
		return t;
	}
	
	// method to compute the differential of f
	private double [] [] Df (double g, double k, double [] v)
	{
		// local variable
		double [] [] vv = new double [v.length] [v.length];
		
		// initialize the components
		vv [0] [0] = g / sq (Math.cosh(g * v[0])) - 1;
		vv [0] [1] = -1.0;
		vv [0] [2] = 0.0;
		
		vv [1] [0] = 1.0; vv[1] [1] = -(k + 1.0); vv[1] [2] = 0.0;
		
		vv [2] [0] = v[2] / sq(Math.cosh(g * v[0]));
		vv [2] [1] = -v[2];
		vv [2] [2] = -v[1] + Math.tanh(g * v[0]);
		
		return vv;
	}
	
	// method to define the function for the Newton-Raphson iteration
	private double [] f (double g, double k, double [] v)
	{
		// local variable
		double [] vv = new double [v.length];
		
		// initailize components
		vv[0] = -v[1] + Math.tanh(g * v[0]) - v[0];
		vv[1] = v[0] - (k + 1) * v[1];
		vv[2] = -(v[1] + Math.tanh(g * v[0])) * v[2];
		
		return vv;
	}
	
	// method to perform square of values
	private double sq (double x) { return x * x; }
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new InverseTest ();
		System.exit(0);
	}

	// INNER CLASS
	private class MatrixException extends Exception
	{	
		public MatrixException (String message) {
			super (message);
		}
	}
}
