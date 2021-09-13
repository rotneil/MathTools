package math.tools;

public final class Matrix {
	
	public static void main (String [] args) {
		// testing the matrix product
		double [] x = Matrix.product(
				new double [] [] {{2, 3, 1}, {4, 5, 1}, {1, 1, 1}},
				new double [] {1, 2, 3});
		for (int i = 0; i < x.length; ++i)
			System.out.println ("x" + (i + 1) + " = " + x[i]);
		
		/*
		// Testing LUSolver
		double [] x = Matrix.LUSolver(
				new double [] [] {{2.0, 3.0, -4.0}, {1.0, 1, 1}, {5, 10.0, -1}}, 
				new double [] {-0.8, 0.75, -1.75});
		
		for (int i = 0; i < x.length; ++i)
			System.out.println ("x" + (i + 1) + " = " + x[i]);
		
		/* Testing the Lyapunove/variation library
		StringBuffer [] f = new StringBuffer [] {
				new StringBuffer ("s * (y - x)"),
				new StringBuffer ("x * (r - z) - y"),
				new StringBuffer ("x * y - b * z")
		};
		
		String [] x = new String [] {"x", "y", "z"};
		String [] var = new String [] {"x", "y", "z", "s", "r", "b"};
		double [] val = new double [] {5.0, 5.0, 5.0, 10.0, 28.0, (8.0 / 3)};
		
		ExpNode [] fx = Matrix.infixToExpNode(f, var);
		
		double [] h = Matrix.getLyapunoveExpCont(fx, x, 200000, var, val, 0.01);
		
		Matrix.display(h);*/
	}
	/**
	 * This method forms the chi function X(u) for the stability condition
	 * from the Jacobian expression Df where the eigen value is given by u. 
	 * This is done by replacing the diagonal elements by a new BinOpNode.
	 * 
	 * @param fx The state function as an ExpNode
	 * @param x The independent coordinate of the system
	 * @param re The real part of the complex eigenvalue
	 * @param im The imaginary part of the complex eigenvalue
	 * @return The chi function X as a DetNode
	 */
	public static DetNode getX (ExpNode [] fx, String [] x, String re, String im)
	{
		return new DetNode (fx, x, re, im);
	}
	
	
	
	/**
	 * This method returns a PeriodicDetNode that can be used in the 
	 * numerical analysis of periodic limit cycle bifurcation computations.
	 * 
	 * @param x The independent coordinates of the state
	 * @param re The real part of the complex eigenvalue
	 * @param im The imaginary part of the complex eigenvalue
	 * @param param The string value of the parameter
	 * @return PeriodicDetNode
	 */
	public static PeriodicDetNode getX (String [] x, String re, String im, 
			String param)
	{
		return new PeriodicDetNode (x, re, im, param);
	}
	
	/**
	 * This method converts an equation expression to its equivalent
	 * ExpNode which can be computed or differentiated. The method first
	 * converts the infix expression to postfix and then to ExpNode.
	 * @param infix The array of equation without the equality sign
	 * @param variables The variables with which the expressions are built
	 * @return
	 */
	public static ExpNode [] infixToExpNode (
			StringBuffer [] infix, String [] variables)
	{
		// instantiate the return
		ExpNode [] ffx = new ExpNode [infix.length];
		
		// convert the equation to postfix
		for (int i = 0; i < infix.length; ++i)
			ffx[i] = Algebra.infixToBinOpNode(infix[i], variables);
		
		return ffx;
	}
	
	/**
	 * This method takes the equation expression and the list of the
	 * independent variables, and returns the first variation equation. The variation equation
	 * are returned as a one dimensional array.
	 * 
	 * If the Jacobian matrix is J(x) and the variational matrix is
	 * 						| x11	x12	...	x1n |
	 * 			A(t) = J(x)	| x21	x22	...	x2n |
	 * 						|  .	 .		.	|
	 * 						| xn1	xn2	...	xnn |
	 * Then the variation is given as
	 * 
	 * 						| x11	x12	...	x1n |
	 * 		Variation (t) =	| x21	x22	...	x2n |
	 * 						|  .	 .		.	|
	 * 						| xn1	xn2	...	xnn |
	 * 	 
	 * 
	 * @param fx The expression node
	 * @param x The independent variables
	 * @return The state equation, and the variation equation. The array is
	 * arranged as {f1(x), f2(x), ..., fn(x), x11, x12, ..., 
	 * x1n, x21, x22, ..., x2n, ..., xn1, xn2, ..., xnn}
	 */
	public static ExpNode [] getVariationEq(ExpNode [] fx, String [] x)
	{
		// get the Jacobian matrix
		int n = x.length;
		ExpNode [] first = new ExpNode [n * (n + 1)];
		ExpNode [] [] Jx = jacobian (fx, x);
		ExpNode [] [] Vx = new ExpNode [n][n];
		ExpNode [] p = new ExpNode [n];		// used for proper variation
		
		// the proper variational equation
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				// continue with the proper variation
				for (int k = 0; k < n; ++k)
					p [k] = BinOpNode.make ('*', Jx[i][k], 
							new VariableNode("x" + (k + 1) + (j + 1)));
				
				Vx [i] [j] = p [0];
				for (int l = 1; l < n; ++l)
					Vx [i][j] = BinOpNode.make ('+', Vx[i][j], p[l]);
			}
		}
		
		// linearize the variation equations
		for (int i = 0; i < n; ++i)
			first [i] = fx[i];
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				first [n * (i + 1) + j] = Vx [i][j];
		
		return first;
	}
	
	/**
	 * This method returns the independent variables of the Variation 
	 * equations returned by the method getVariationEquation (fx, x).
	 * 
	 * @param x
	 * @return
	 */
	public static String [] getVariationIndependentVar (String [] x)
	{
		int n = x.length;
		String [] var = new String [n * (n + 1)];
		
		// put the x variables
		for (int i = 0; i < n; ++i)
			var[i] = x[i];
		
		// put those of the proper variation equations
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				var [n * (i + 1) + j] = "x" + (i + 1) + (j + 1);
		
		return var;
	}
	
	/**
	 * This method returns the independent variables of the Variation 
	 * equations returned by the method getFirstVariationEquation (fx, x).
	 * 
	 * @param var The original model variables
	 * @return The combination of the model variables with those of the 
	 * variation equation
	 */
	public static String [] getVariationVar (String [] var, String [] x)
	{
		int n = x.length;
		String [] mVar = new String [n * n + var.length];
		
		// put the element variables
		for (int i = 0; i < var.length; ++i)
			mVar[i] = var[i];
		
		// put those of the proper variation equations
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				mVar [var.length + n * i + j] = "x" + (i + 1) + (j + 1);
		
		return mVar;
	}
	
	/**
	 * This method is used to initialize the values for the elements of the
	 * first variation equation which does not include the parameter variation. 
	 * The double array values val are included in the returned array.
	 *
	 * @param val The initial values without those of the variation equation
	 * @param x The independent variables of the system.
	 * @return The double arrays with the elements of the variation equation
	 */
	public static double [] getVariationInitVal (double [] val, String [] x)
	{
		int n = x.length;
		double [] init = new double [val.length + n * n];
		
		// copy contents of the val into init
		for (int i = 0; i < val.length; ++i) init [i] = val [i];
		
		// copy contents of identity into init
		for (int i = 0; i < n; ++i)
			init [val.length + i * (n + 1)] = 1.0;
		
		return init;
	}
	
	
	/**
	 * This method takes the equation expression and the list of the
	 * independent variables, and returns the first variation equation. This
	 * equation includes the equation for the state. The variation equation
	 * are returned as a one dimensional array.
	 * 
	 * If the Jacobian matrix is J(x) and the variational matrix is
	 * 						| x11	x12	...	x1n |
	 * 			A(t) = J(x)	| x21	x22	...	x2n |
	 * 						|  .	 .		.	|
	 * 						| xn1	xn2	...	xnn |
	 * Then the variation is given as
	 * 
	 * 						| x11	x12	...	x1n |
	 * 		Variation (t) =	| x21	x22	...	x2n |
	 * 						|  .	 .		.	|
	 * 						| xn1	xn2	...	xnn |
	 * 	 
	 * 
	 * @param fx The expression node
	 * @param x The independent variables
	 * @param param
	 * 
	 * @return The state equation, and the variation equation. The array is
	 * arranged as {f1(x), f2(x), ..., fn(x), x11, x12, ..., x1n, x21, x22, 
	 * ..., x2n, ..., xn1, xn2, ..., xnn, kk1, kk2}
	 */
	public static ExpNode [] getFirstVariationEq(ExpNode [] fx, 
			String [] x, String param)
	{
		// get the Jacobian matrix
		int n = x.length;
		ExpNode [] first = new ExpNode [n * (n + 2)];
		ExpNode [] [] Jx = jacobian (fx, x);
		ExpNode [] [] Vx = new ExpNode [n][n];
		ExpNode [] Kx = new ExpNode [n];
		ExpNode [] p = new ExpNode [n];		// used for proper variation
		ExpNode [] q = new ExpNode [n];		// used for parameter variation
		
		// the proper variational equation
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				// the parameter variational equation
				q [j] = BinOpNode.make ('*', Jx[i][j], 
							new VariableNode ("x" + (j + 1) + "k"));
				
				// continue with the proper variation
				for (int k = 0; k < n; ++k)
					p [k] = BinOpNode.make ('*', Jx[i][k], 
							new VariableNode("x" + (k + 1) + (j + 1)));
				
				Vx [i] [j] = p [0];
				for (int l = 1; l < n; ++l)
					Vx [i][j] = BinOpNode.make ('+', Vx[i][j], p[l]);
			}
			Kx [i] = fx[i].derivative(param);
			for (int k = 0; k < n; ++k)
				Kx[i] = BinOpNode.make ('+', q[k], Kx[i]);
		}
		
		// linearize the variation equations
		for (int i = 0; i < n; ++i)
			first [i] = fx[i];
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				first [n * (i + 1) + j] = Vx [i][j];
		for (int i = 0; i < n; ++i)
			first[n * (n + 1) + i] = Kx[i];
		
		return first;
	}
	
	/**
	 * This method returns the independent variables of the Variation 
	 * equations returned by the method getFirstVariationEquation (fx, x).
	 * 
	 * @param x
	 * @return
	 */
	public static String [] getFirstVariationIndependentVar (String [] x)
	{
		int n = x.length;
		String [] var = new String [n * (n + 2)];
		
		// put the x variables
		for (int i = 0; i < n; ++i)
			var[i] = x[i];
		
		// put those of the proper variation equations
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				var [n * (i + 1) + j] = "x" + (i + 1) + (j + 1);
		
		// put those of the parameter variation
		for (int i = 0; i < n; ++i)
			var [n * (n + 1) + i] = "x" + (i + 1) + "k";
		
		return var;
	}
	
	/**
	 * This method returns the independent variables of the Variation 
	 * equations returned by the method getFirstVariationEquation (fx, x).
	 * 
	 * @param var The original model variables
	 * @param x The independent variables
	 * @return The combination of the model variables with those of the 
	 * variation equation
	 */
	public static String [] getFirstVariationVar (String [] var, String [] x)
	{
		int n = x.length;
		String [] mVar = new String [n * (n + 1) + var.length];
		
		// put the element variables
		for (int i = 0; i < var.length; ++i)
			mVar[i] = var[i];
		
		// put those of the proper variation equations
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				mVar [var.length + n * i + j] = "x" + (i + 1) + (j + 1);
		
		// put those of the parameter variation
		for (int i = 0; i < n; ++i)
			mVar [var.length + n * n + i] = "x" + (i + 1) + "k";
		
		return mVar;
	}
	
	/**
	 * This method is used to initialize the values for the elements of the
	 * first variation equation. The double array values val are included in 
	 * the returned array.
	 * 
	 * This method only initializes the indices of the diagonal elements of 
	 * identity matrix.
	 * 
	 * @param val The initial values without those of the variation equation
	 * @param x The independent variables of the system.
	 * @return The double arrays with the elements of the variation equation
	 */
	public static double [] getFirstVariationInitVal (double [] val, String [] x)
	{
		int n = x.length;
		double [] init = new double [val.length + n * (n + 1)];
		
		// copy contents of the val into init
		for (int i = 0; i < val.length; ++i) init [i] = val [i];
		
		// copy contents of identity into init
		for (int i = 0; i < n; ++i)
			init [val.length + i * (n + 1)] = 1.0;
		
		return init;
	}
	
	/**
	 * This method returns the matrix form of the first variation equation
	 * independent variables as string
	 * 
	 * @param x The independent variables
	 * @return
	 */
	public static ExpNode [] [] getFirstVariationXNode (String [] x)
	{
		int n = x.length;
		ExpNode [] [] wx = new ExpNode [n] [n];
		
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				wx [i][j] = ExpNode.makeVariableNode("x" + (i + 1) + (j + 1));
		
		return wx;
	}
	
	/**
	 * This method returns the variables of the first variation parameter
	 * coordinates as a matrix of order n X 1.
	 * 
	 * @param x
	 * @return
	 */
	public static ExpNode [] getFirstVariationParamNode (String [] x)
	{
		int n = x.length;
		ExpNode[] k = new ExpNode [n];
		
		for (int i = 0; i < n; ++i)
			k[i] = new VariableNode ("x" + (i + 1) + "k");
		
		return k;
	}
	
	/**
	 * This method returns the equations required for the second
	 * variational equation computation. The equations are returned as a
	 * one-dimensional array.
	 * 
	 * If the first variational matrix is given as
	 * 		| dw1/dx01		dw1/dx02	...		dw1/dx0n	|	  | dw1/dk |
	 * 		| dw2/dx01		dw2/dx02	...		dw2/dx0n	|	  | dw2/dk |
	 * 		| 	...			   ...		...			...		| and |   ...  |
	 * 		| dwn/dx01		dwn/dx02	...		dwn/dx0n	|	  | dwn/dk |
	 * 
	 * and is noted as
	 * 		| x11   x12   x13   ...   x1n |		| x1k |
	 * 		| x21   x22   x23   ...   x2n |		| x2k |
	 * 		| ...	...   ...	...	  ... | and | ... |
	 * 		| xn1   xn2   xn3   ...   xnn |		| xnk |
	 * 		
	 * while the second variational matrix is given as
	 * 		| d2w1/dx01^2	d2w1/dx01dx02	d2w1/dx01dx03	... 
	 * 			d2w1/dx01dx0n	d2w1/dx02^2		d2w1/dx02dx03	...	d2w1/dx0n^2	|
	 * 		| d2w2/dx01^2	d2w2/dx01dx02	d2w2/dx01dx03	... 
	 * 			d2w2/dx01dx0n	d2w2/dx02^2		d2w2/dx02dx03	...	d2w2/dx0n^2	|
	 * 		|	...				...				...			...					|
	 * 		| d2wn/dx01^2	d2wn/dx01dx02	d2wn/dx01dx03	... 
	 * 			d2wn/dx01dx0n	d2wn/dx02^2		d2wn/dx02dx03	...	d2wn/dx0n^2|
	 * 
	 * and they are denoted as
	 * 
	 * 	| x1x11  x1x12  x1x13 ... x1x1n ... x1x22 x1x23 ... x1x2n ... x1xnn |
	 * 	| x2x11  x2x12  x2x13 ... x2x1n ... x2x22 x2x23 ... x2x2n ... x2xnn |
	 * 	|  ...    ...    ...       ...  ...  ...   ...       ...	   ...  |
	 * 	| x1x11  x1x12  x1x13 ... x1x1n ... x1x22 x1x23 ... x1x2n ... x1xnn |
	 * 
	 * and that of the parameter as
	 * 
	 * 	| d2w1/dx01dk   d2w1/dx02dk ... d2w1/dx0ndk |     | x11k x12k ... x1nk |
	 *  | d2w2/dx01dk   d2w2/dx02dk ... d2w2/dx0ndk |     | x21k x22k ... x2nk |
	 *  |	 ...		   ...				...     | as  |  ...  ...	   ... |
	 *  | d2wn/dx01dk   d2wn/dx02dk ... d2wn/dx0ndk |     | xn1k xn2k ... xnnk |
	 * 
	 * Then the return linear equation is in this order
	 *  {x1  x2 ... xn  x11  x12 ... x1n  x21  x22 ... x2n ... xn1  xn2  ...
	 *  	xnn  x1k  x2k ... xnk  x1x11  x1x12  x1x13 ... x1x1n  x1x22  
	 *  	x1x23 ... x1x2n ... x1xnn  x2x11  x2x12  x2x13 ... x2x1n x2x22  
	 *  	x2x23 ... x2x2n ... x2xnn ... xnx11  xnx12  xnx13 ... xnx1n  xnx22
	 *  	xnx23 ... xnx2n ... xnxnn  x11k  x12k  x13k ... x1nk  x21k  x22k
	 *  	x23k  ... x2nk  ... xn1k   xn2k  xn3k ... xnnk}
	 *  
	 * @param fx The state equations
	 * @param x The independent variables of the state equation
	 * @param param The parameter required for the second variation
	 * 
	 * @return
	 */
	public static ExpNode [] getSecondVariationEq (ExpNode [] fx,
			String [] x, String param)
	{
		// local variable that has all the terms of the variational equation
		int n = x.length;
		ExpNode [] Vx = new ExpNode [(n * (n + 1) * (n + 4) / 2)];
		
		// copy the content of the first variation into Vx
		ExpNode [] first = getFirstVariationEq (fx, x, param);
		for (int i = 0; i < first.length; ++i)
			Vx [i] = first [i];
		
		// begin the evaluation of the second variation equation
		// by calling for the jacobian of fx
		ExpNode [] [] jx = jacobian (fx, x);
		ExpNode [] [] v2x = new ExpNode [n][(n * (n + 1) / 2)];
		ExpNode [] [] v2k = new ExpNode [n][n];
		ExpNode sum;
		ExpNode [] jxSum = new ExpNode [n];		// the sum of product
		ExpNode [] j2x, jfk;			// the product
		
		int index;
		for (int i = 0; i < n; ++i) {
			index = 0;
			for (int j = 0; j < n; ++j) {
				for (int k = j; k < n; ++k) {
					// perform the first part of the second variational eq
					// that is, (Dfx)i * {d2wi/dx0jdx0k}
					sum = BinOpNode.make ('*', jx [i][0], 
							new VariableNode ("x1" + 
									"x" + (j + 1) + (k + 1)));
					for (int l = 1; l < n; ++l)
						sum = BinOpNode.make ('+', sum, BinOpNode.make ('*', 
							jx[i][l], new VariableNode ("x" + (l + 1) + 
									"x" + (j + 1) + (k + 1))));
					
					// now take each element of the jacobian and differentiate
					// with respect to the independent variables and multiply
					// with appropriate element of the first variational eq
					for (int l = 0; l < n; ++l) {
						j2x = derivate (jx[i][l], x);
						
						// post multiply the elements of the column of first
						// variational equation
						for (int m = 0; m < n; ++m)
							j2x[m] = BinOpNode.make ('*', j2x[m], 
								new VariableNode ("x" + (m + 1) + (j + 1)));
								
						// perform the sum of the product of fx.der * x11
						jxSum [l] = j2x[0];
						for (int m = 1; m < n; ++m)
							jxSum[l] = BinOpNode.make ('+', jxSum[l], j2x[m]);
					}
					
					// now form the sum of Dfx * {d2wk/dx0j}
					for (int l = 0; l < n; ++l)
						sum = BinOpNode.make ('+', sum, BinOpNode.make ('*', 
							jxSum[l], new VariableNode (
								"x" + (l + 1) + (k + 1))));
					
					// assign sum into the second variational matrix
					v2x[i][index++] = sum;
				}
			}
		}
		
		// perform the second variational equation for the parameters
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				// start Dfx * { d2wi/dxj_dk }
				sum = BinOpNode.make ('*', jx[i][0], new VariableNode (
						"x1" + (j + 1) + "k"));
				for (int k = 1; k < n; ++k)
					sum = BinOpNode.make ('+', sum, BinOpNode.make ('*', jx[i][k], 
						new VariableNode ("x" + (k + 1) + (j + 1) + "k")));
				
				// start d/dx0j { Dfx } * { dwi/dk }
				for (int k = 0; k < n; ++k) {
					// derivate each element of jx wrt to x
					j2x = derivate (jx[i][k], x);
					
					// pre-multiply contents of j2x with xji
					for (int l = 0; l < n; ++l)
						j2x[l] = BinOpNode.make ('*', j2x[l], 
							new VariableNode ("x" + (l + 1) + (j + 1)));
					
					// perform the sum
					jxSum[k] = j2x[0];
					for (int l = 1; l < n; ++l)
						jxSum[k] = BinOpNode.make ('+', jxSum[k], j2x[l]);
					
					// finish the operation d/dx0i (Dfx)j * {dwi/dk}
					jxSum[k] = BinOpNode.make ('*', jxSum[k],
						new VariableNode ("x" + (k + 1) + "k"));
				}
				
				// start d/dx0i { df } = d2f/dxdk * { dw/dxi }
				jfk = derivate (fx[i].derivative(param), x);
				
				// implement the product
				for (int k = 0; k < jfk.length; ++k)
					jfk [k] = BinOpNode.make ('*', jfk [k], 
							new VariableNode ("x" + (k + 1) + (j + 1)));
				
				// now perform the row sum
				for (int k = 0; k < n; ++k)
					sum = BinOpNode.make ('+', sum, BinOpNode.make ('+',
							jxSum[k], jfk [k]));
				
				v2k [i][j] = sum;
			}
		}
		
		// linearize the elements of v2x into Vx
		int xStart = n * (n + 2);
		int xlen = n * (n + 1) / 2;
		int kStart = n * (2 * n + 4 + n * (n + 1)) / 2;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < v2x[i].length; ++j)
				Vx[xStart + i * xlen + j] = v2x[i][j];
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				Vx[kStart + i * n + j] = v2k[i][j];
		
		// return the second variation equation
		return Vx;
	}
	
	/**
	 * @see {@link #getSecondVariationEq(ExpNode[], String[], String)}
	 * @param x
	 * @return
	 */
	public static String [] getSecondVariationIndependentVar (String [] x)
	{
		int n = x.length;
		String [] variable = new String [(n * (n + 1) * (n + 4) / 2)];
		
		// copy terms of the first variation equation
		String [] first = getFirstVariationIndependentVar (x);
		for (int i = 0; i < first.length; ++i)
			variable [i] = first [i];
		
		int start = n * (n + 2);
		int len = n * (n + 1) / 2;
		int index;
		// copy the second term
		for (int i = 0; i < n; ++i) {
			index = 0;
			for (int j = 0; j < n; ++j) {
				for (int k = j; k < n; ++k) 
					variable [start + i * len + index++] = "x" + (i + 1) +
						"x" + (j + 1) + (k + 1);
			}
		}
		
		// copy terms of the parameters
		start += n * n * (n + 1) / 2;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				variable [start + i * n + j] = "x" + (i + 1) + (j + 1) + "k";
		
		return variable;
	}
	
	/**
	 * @see {@link #getSecondVariationEq(ExpNode[], String[], String)}
	 * @param var
	 * @param x
	 * @return
	 */
	public static String [] getSecondVariationVar (String [] var, String [] x)
	{
		int n = x.length;
		String [] first = getFirstVariationVar (var, x);
		String [] variable = new String [first.length + n * n * (n + 3) / 2];
		
		for (int i = 0; i < first.length; ++i)
			variable[i] = first[i];
		
		int start = first.length;
		int len = n * (n + 1) / 2;
		int index;
		// copy the second term
		for (int i = 0; i < n; ++i) {
			index = 0;
			for (int j = 0; j < n; ++j) {
				for (int k = j; k < n; ++k) 
					variable [start + i * len + index++] = "x" + (i + 1) +
						"x" + (j + 1) + (k + 1);
			}
		}
		
		// copy terms of the parameters
		start += n * n * (n + 1) / 2;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				variable [start + i * n + j] = "x" + (i + 1) + (j + 1) + "k";
		
		return variable;
	}
	
	/**
	 * This method is used to initialize the values for the elements of the
	 * second variation equation. The double array values val are included in 
	 * the returned array.
	 *
	 * @see {@link #getFirstVariationInitValues (double[], String[])}
	 * 
	 * @param val The initial values without those of the variation equation
	 * @param x The independent variables of the system.
	 * @return The double arrays with the elements of the variation equation
	 */
	public static double [] getSecondVariationInitVal (double [] val, String [] x)
	{
		int n = x.length;
		double [] first = getFirstVariationInitVal (val, x);
		double [] init = new double [first.length + n * n * (n + 3) / 2];
		
		// copy contents of the val into init
		for (int i = 0; i < first.length; ++i) init [i] = first [i];
		
		return init;
	}
	
	
	/**
	 * This method returns the values of an ExpNode in array
	 * @param fx The ExpNode Array
	 * @param var The variable definition
	 * @param val The corresponding double values
	 * @return
	 */
	public static double [] values (ExpNode [] fx, String [] var, double [] val)
	{
		// local variable
		double [] v = new double [fx.length];
		for (int i = 0; i < fx.length; ++i)
			v [i] = fx[i].value(var, val);
		return v;
	}
	
	/**
	 * This method computes the double value for a 2-dimensional ExpNode array
	 * @param fx The ExpNode array
	 * @param var The variables
	 * @param val variables corresponding values
	 * @return
	 */
	public static double [] [] values (ExpNode[][] fx, 
			String[] var, double[] val)
	{
		// local variable
		double [] [] v = new double [fx.length] [fx[0].length];
		for (int i = 0; i < fx.length; ++i)
			for (int j = 0; j < fx[0].length; ++j)
				v [i][j] = fx[i][j].value(var, val);
		return v;
	}
	
	/**
	 * This method returns the double values of an array of ExpNode fx which 
	 * has an index indicating a real number and an imaginary number. It does
	 * so by identifying which of the equations is a BinOpNode or a DetNode, 
	 * and then evaluates its ComplexNumber value in order to appropriately
	 * return the Real Index and the Imaginary Index.
	 * 
	 * @param fx ExpNode equation array
	 * @param var The variables of the ExpNode
	 * @param val The corresponding values
	 * @param realIndex The index of the Real Chi function
	 * @param imaginaryIndex The imaginary of the Chi function
	 * @return
	 * @throws MatrixException
	 */
	public static double [] values (ExpNode [] fx, String [] var, double [] val,
			int realIndex, int imaginaryIndex) throws MatrixException
	{
		// check that indices for real and imaginary equations are valid
		if (realIndex < 0 || realIndex >= fx.length)
			throw new MatrixException ("Invalid index for the real " +
					"equation. Index = " + realIndex);
		if (imaginaryIndex < 0 || imaginaryIndex >= fx.length)
			throw new MatrixException ("Invalid index for the " +
					"imaginary equation. Index = " + imaginaryIndex);
		
		double [] v = new double [fx.length];
		
		// evaluate the non complex functions
		for (int i = 0; i < realIndex; ++i)
			v[i] = fx[i].value(var, val);
		
		// Express the complex index as X
		ExpNode X = fx [realIndex];
		ComplexNumber complex = X.complexValue(var, val);
		
		// instantiate the returned double array
		v[realIndex] = complex.getRealNumber();
		v[imaginaryIndex] = complex.getImaginaryNumber();
		return v;
	}
	
	/**
	 * This method returns the double values of an array of ExpNode fx which 
	 * has an index indicating a real number and an imaginary number. It does
	 * so by identifying which of the equations is a BinOpNode or a DetNode, 
	 * and then evaluates its ComplexNumber value in order to appropriately
	 * return the Real Index and the Imaginary Index.
	 * 
	 * @param fx ExpNode equation array
	 * @param var The variables of the ExpNode
	 * @param val The corresponding values
	 * @param realIndex The index of the Real Chi function
	 * @param imaginaryIndex The imaginary of the Chi function
	 * @return
	 * @throws MatrixException
	 */
	public static double [] [] values (ExpNode[][] fx, String [] var,
		double [] val, int realIndex, int imaginaryIndex) throws MatrixException
	{
		// check that indices for real and imaginary equations are valid
		if (realIndex < 0 || realIndex >= fx.length)
			throw new MatrixException ("Invalid index for the real " +
					"equation. Index = " + realIndex);
		if (imaginaryIndex < 0 || imaginaryIndex >= fx.length)
			throw new MatrixException ("Invalid index for the " +
					"imaginary equation. Index = " + imaginaryIndex);
		
		// define the return value
		double [] [] v = new double [fx.length][fx[0].length];
		
		// evaluate values for the nonComplex equations
		for (int i = 0; i < realIndex; ++i)
			v[i] = values (fx[i], var, val);
		
		// compute values for the complex equations by taking one of 
		// complex chi equations and putting the result for the cases
		// of the real and the imaginary
		ExpNode [] X = fx [realIndex];
		ComplexNumber complex;
		for (int i = 0; i < X.length; ++i) {
			complex = X[i].complexValue(var, val);
			
			v[realIndex][i] = complex.getRealNumber();
			v[imaginaryIndex][i] = complex.getImaginaryNumber();
		}
		
		return v;
	}
	
	/**
	 * This method evaluates the Jacobian matrix of an array of ExpNode using
	 * the array of variables x as the independent variables
	 * @param fx The array of the dependent variables in ExpNode
	 * @param x The array of independent variables
	 * @return
	 */
	public static ExpNode [] [] jacobian (ExpNode [] fx, String [] x)
	{
		// instantiate the jacobian matrix
		ExpNode [] [] dfx = new ExpNode [fx.length] [x.length];
		
		// differentiate each function fx wrt to x
		for (int i = 0; i < fx.length; ++i)
			dfx[i] = derivate(fx[i], x);
		
		return dfx;
	}
	
	/**
	 * This method returns the derivative of equation fx wrt to each
	 * of the element of the x. That is, {df(x1)  df(x2) ... df(xn) }
	 * 
	 * @param fx The equation
	 * @param x The independent variables
	 * @return
	 */
	public static ExpNode [] derivate (ExpNode fx, String [] x)
	{
		ExpNode [] dfx = new ExpNode [x.length];
		
		for (int i = 0; i < x.length; ++i)
			dfx[i] = fx.derivative(x[i]);
		
		return dfx;
	}
	
	/**
	 * This method returns the eigen values of Matrix a.
	 * @see Polynomial.getRoots(a, n);
	 * @see Matrix.getCharacteristicsEquation (a)
	 * @param a
	 * @return
	 * @throws MatrixException
	 */
	public static ComplexNumber [] getEigenValues (double [] [] a)
		throws MatrixException
	{
		return Algebra.getRoots(Matrix.getCharacteristicEquation(a),
				a.length);
	}
	
	/**
	 * This method uses Faddeev-Leverrier's method of Characteristic equation
	 * evaluation to derive the characteristics polynomial from a square matrix.
	 * The returned array has the coefficients of the matrix as
	 * 		f(x) = a0 + a1 x + a2 + x^2 + ... + an+1 x^(n+1) + an x^n. 
	 * @param a
	 * @return
	 * @throws MatrixException
	 */
	public static double [] getCharacteristicEquation (double [] [] a)
		throws MatrixException
	{
		// local variables
		int n = a.length;
		double [] p = new double [n];
		double [] an = new double [n + 1];
		
		// assign the initial values of b and p
		double [] [] b = a.clone();
		p[0] = Matrix.getTrace (b);
		
		// assign subsequent iterate
		for (int k = 1; k < n; ++k) {
			b = Matrix.product(a, Matrix.subtract(b, 
					Matrix.product(p[k - 1], Matrix.getIdentity(a.length))));
			// assign the new p
			p[k] = Matrix.getTrace(b) / (k + 1);
		}
		
		// interchange the coefficients of p
		for (int i = n; i > 0; --i)
			an[n - i] = -p[i - 1];
		an[n] = 1.0;
		
		return an;
	}
	
	/**
	 * This method returns the trace of a matrix A which is the sum of the 
	 * diagonal elements of matrix A
	 * @param a
	 * @return
	 */
	public static double getTrace (double [] [] a) 
	{
		double tr = 0.0;
		for (int i = 0; i < a.length; ++i)
			tr += a[i] [i];
		return tr;
	}
	
	/**
	 * This method uses Gauss elimination method to compute the determinant of a
	 * matrix a.
	 * @param a
	 * @return
	 */
	public static ComplexNumber getDeterminant (ComplexNumber [] [] a) 
			throws MatrixException
	{
		// check for zero rows
		for (int row = 0; row < a.length; ++row)
			if (checkForZeroRow (a[row]))
				return new ComplexNumber ();
		
		// local variables
		int [] pivot = new int [a.length];
		double tol = 1.0E-12;
		int [] multiplier = new int [a.length];
		for (int i = 0; i < a.length; ++i) multiplier [i] = 1;
		
		// decompose matrix a into lower and upper triangular matrices
		ComplexNumber [] [] d = decompose (a, pivot, 
				new ComplexNumber [a.length], tol, multiplier);
		
		// compute the determinant as the product of the diagonal elements
		ComplexNumber det = d[pivot[0]][0].product(multiplier[0]);
		for (int i = 1; i < d.length; ++i)
			det = det.product(d[pivot[i]][i].product(multiplier[i]));
		
		return det;
	}
	
	/**
	 * This method checks the row if all the elements are zero.
	 * @param row
	 * @return True if it has rows which are all zero
	 */
	public static boolean checkForZeroRow (ComplexNumber [] row)
	{
		for (int col = 0; col < row.length; ++col)
			if (row[col].getRealNumber() != 0.0 || 
				row[col].getImaginaryNumber() != 0.0)
				return false;
		return true;
	}
	
	/**
	 * This method uses Gauss elimination method to compute the determinant of a
	 * matrix a.
	 * @param a
	 * @return
	 */
	public static double getDeterminant (double [] [] a) throws MatrixException
	{
		// check for zero rows
		for (int row = 0; row < a.length; ++row)
			if (checkForZeroRow (a[row]))
				return 0.0;
		
		// local variables
		int [] pivot = new int [a.length];
		double tol = 1.0E-12;
		int [] multiplier = new int [a.length];
		for (int i = 0; i < a.length; ++i) multiplier [i] = 1;
		
		// decompose matrix a into lower and upper triangular matrices
		double [] [] d = decompose (a, pivot, new double [a.length], tol, multiplier);
		
		// compute the determinant as the product of the diagonal elements
		double det = 1.0;
		for (int i = 0; i < d.length; ++i)
			det *= d[pivot[i]][i] * multiplier[i];
		
		return det;
	}
	
	/**
	 * This method checks the row if all the elements are zero.
	 * @param row
	 * @return True if it has rows which are all zero
	 */
	public static boolean checkForZeroRow (double [] row)
	{
		for (int col = 0; col < row.length; ++col)
			if (row[col] != 0.0)
				return false;
		return true;
	}
	
	/**
	 * This checks the rows of either a DetNode or a PeriodicNode for instance
	 * of a zero-element row.
	 * 
	 * @param X The DetNode or PeriodicDetNode.
	 * @return True, if any row has zero-element.
	 */
	public static boolean hasZeroRow (ExpNode [][] X)
	{
		boolean zeroRow = true;
		// check rows
		for (int i = 0; i < X.length; ++i) {
			zeroRow = true;
			
			// loop throw the columns of the rows
			for (int j = 0; j < X[i].length; ++j) {
				if (X[i][j] instanceof ConstantNode) {
					if (((ConstantNode) X[i][j]).number != 0.0) {
						zeroRow = false;
						break;
					}
				}
				else {
					zeroRow = false;
					break;
				}
			}
			if (zeroRow)
				return zeroRow;
		}
		return zeroRow;
	}
	
	/**
	 * This method uses Gauss elimination method to seperate matrix A into its corresponding
	 * lower and upper triangular matrices. The elements of the formed upper triangular matrix
	 * represent the U(i,j) from i = 0 to n and j = i to n, while the lower triangular matrix
	 * are L(i, j) from i = 1 to n for j = 0 to i.
	 * Note: The diagonal element of L are all 1's- so there's no need for representing them
	 * in a real matrix.
	 * 
	 * @param a The matrix to be decomposed into upper and lower matrices
	 * @param pivot The pivot transformation order. During decomposition, rows are supposed to
	 * be interchanged to enhance the elimination process. But instead, rows are retained but 
	 * the indices of the pivot are stored in this argument
	 * @param scale The scaling factor of the matrix. This is the maximum element in a row
	 * @param tol The tolerance, below which an element is assumed to be zero.
	 * @param multiplier Integer that tracks row interchange for computation of 
	 * the determinant of a matrix
	 * @return The decomposed matrix that has the upper and lower triangular matrix
	 * @throws MatrixException
	 */
	public static ComplexNumber [] [] decompose (ComplexNumber [] [] a, 
			int [] pivot, ComplexNumber [] scale, double tol, 
			int [] multiplier)	throws MatrixException
	{
		// confirm the validity of matrix a
		if (a.length != a[0].length)
			throw new MatrixException ("Matrix a is not a square matrix");
		if (pivot.length != a.length)
			throw new MatrixException ("The rows of matrix a must equal to the length of pivot");
		if (scale.length != a.length)
			throw new MatrixException ("The rows of matrix a must equal to the length of scale");
		
		// local variable
		ComplexNumber [] [] c = new ComplexNumber [a.length][a.length];
		for (int i = 0; i < a.length; ++i)
			for (int j = 0; j < a[i].length; ++ j)
				c[i][j] = a[i][j];
		
		// set the scale and the pivot to reflect the normal indices
		setScale (c, pivot, scale);
		
		// now begin the elimination process
		for (int k = 0; k < c.length - 1; ++k) {
			// call method setPivot
			setPivot (c, pivot, scale, k, multiplier);
			
			// check the system's tolerance before elimination
			ComplexNumber scaledValue = c[pivot[k]][k].divide(scale[pivot[k]]);
			if (Math.abs(scaledValue.getModulus()) < tol)
				throw new MatrixException ("Scaled value " + scaledValue + 
						" is almost equal to zero: Less than the required tolerance " + tol);
			
			// operate on the lower rows
			for (int i = k + 1; i < c.length; ++i) {
				// evaluate the multiplication factor
				ComplexNumber factor = c [pivot[i]][k].divide(c [pivot[k]][k]);
				
				// set the lower triangular matrix element
				c[pivot[i]][k] = factor;
				
				// eliminate pivot row from the current row
				for (int j = k + 1; j < c.length; ++j)
					c[pivot[i]][j] = c [pivot[i]][j].subtract(
							factor.product(c [pivot[k]][j]));
			}
		}
		
		return c;
	}
	
	/**
	 * This method sets the scale of the row operation during Gaussian elimination
	 * method. It seeks in each rows the column that has the highest element.
	 * @param c The matrix whose scales are to be set
	 * @param pivot	The pivot of the matrix c
	 * @param scale The values in the scale
	 */
	private static void setScale (ComplexNumber [] [] c, int [] pivot, 
			ComplexNumber [] scale)
	{
		// iterate through the rows
		for (int i = 0; i < c.length; ++i) {
			pivot[i] = i;
			
			// set the first element of the row to scale
			scale [i] = c[i][0];
			
			// choose the maximum element in the row as the scale
			for (int j = 1; j < c.length; ++j)
				if (c[i][j].compareTo(scale [i]) == 1)
					scale [i] = c[i][j];
		}
	}
	
	/**
	 * This method sets the pivot for the Gauss elimination method.
	 * @param a This is the Matrix to be be decomposed by Gauss elimination method
	 * @param pivot This is the row transformation sequence that is to be set for this row operation
	 * @param scale The scaling value for the rows
	 * @param index The current Row operation that requires pivoting
	 * @param multiplier Integer that tracks row interchanger, used in computation of 
	 * the determinant
	 */
	private static void setPivot (ComplexNumber [] [] a, int [] pivot, 
			ComplexNumber [] scale,	int index, int [] multiplier)
	{
		// check if index is the last row of the array
		if ((index + 1) >= a.length)
			return;
		
		// local variable
		ComplexNumber dummy;
		
		// set the initial value for pivotRow and set the max element
		int rowPivot = index;
 		ComplexNumber max = a[pivot[index]][index].divide(scale[pivot[index]]);
		
		// iterate through the elements of this column for the biggest
		for (int i = index + 1; i < a.length; ++i) {
			// evaluate scaled value of the elements as dummy and compare it with other elements
			dummy = a[pivot[i]][index].divide(scale[pivot[i]]);
			if (dummy.compareTo(max) == 1) {
				max = dummy;
				rowPivot = i; 
			}
		}
		
		// check if rows have not been interchanged
		if (rowPivot == index)
			return;
		
		// interchange the content of pivot at index and rowPivot
		int dummyIndex = pivot [rowPivot];
		pivot [rowPivot] = pivot[index];
		pivot [index] = dummyIndex;
		multiplier [index] *= -1;
	}	// end method setPivot
	
	/**
	 * This method implements the Lower-Upper Triangular matrix decomposition 
	 * to solve a linear system of algebra 
	 * 		[A]{x} = {B}
	 *  
	 * @param A The coefficients of the linear equation of order n X n
	 * @param B The constant of the equation of order n X 1
	 * @return x The value of the column matrix x
	 * 
	 * @throws MatrixException
	 */
	public static ComplexNumber [] LUSolver (ComplexNumber [] [] A, 
			ComplexNumber [] B) throws MatrixException
	{
		// local variable that separates the argument from the return
		ComplexNumber [] x = new ComplexNumber [A[0].length];
		
		// local variables
		ComplexNumber [] scale = new ComplexNumber [A.length];
		int pivot [] = new int [A.length];
		double tol = 1E-18;
		
		// decompose matrix a
		ComplexNumber [] [] lu = decompose(
				A, pivot, scale, tol, new int [A.length]);
		
		// perform the forward and backward substition processes
		substitute(lu, pivot, B, x);
		
		// return the row matrix
		return x;
	}
	
	/**
	 * This method implements the Lower-Upper Triangular matrix decomposition 
	 * to solve a linear system of algebra 
	 * 		[A]{x} = {B}
	 *  
	 * @param A The coefficients of the linear equation of order n X n
	 * @param B The constant of the equation of order n X 1
	 * @return x The value of the column matrix x
	 * 
	 * @throws MatrixException
	 */
	public static double [] LUSolver (double [] [] A, double [] B) 
			throws MatrixException
	{
		// local variable that separates the argument from the return
		double [] x = new double [A[0].length];
		
		// local variables
		double [] scale = new double [A.length];
		int pivot [] = new int [A.length];
		double tol = 1E-18;
		
		// decompose matrix a
		double [] [] lu = decompose(A, pivot, scale, tol, new int [A.length]);
		
		// perform the forward and backward substition processes
		substitute(lu, pivot, B, x);
		
		// return the row matrix
		return x;
	}
	
	/**
	 * This method returns the inverse of matrix Ai such that 
	 * 		Ai x A = In
	 * where In is the identity matrix. This method calls the standard
	 * inverse method {@link(inverse)}
	 * 
	 * @param a
	 * @return
	 * @throws MatrixException
	 */
	public static ComplexNumber [] [] inverse (ComplexNumber [] [] a) 
			throws MatrixException
	{
		return inverse (a, new int [a.length], 
				new ComplexNumber [a.length], 1.0E-18);
	}
	
	
	/**
	 * This method uses Gauss Elimination method to determine the inverse of a 
	 * square matrix [A].
	 * @param A The matrix whose inverse is desired
	 * @param pivot The transformation sequence of the rows during pivot operation
	 * @param scale The scaling factor to enhance elimination process
	 * @param tol The tolerance, which should be as low as possible to avoid division by zero
	 * @return The inverse of matrix A i.e [A]'
	 * @throws MatrixException If any error occurs, especially if a value less than tol is obtained
	 */
	public static ComplexNumber [] [] inverse (ComplexNumber [] [] a, int [] pivot, 
			ComplexNumber [] scale, double tol)	throws MatrixException
	{
		// local variable that separates the argument from the return
		ComplexNumber [] [] x = new ComplexNumber [a.length][a.length];
		ComplexNumber [] b = new ComplexNumber [a.length];
		
		// decompose matrix a
		ComplexNumber [] [] lu = decompose(
				a, pivot, scale, tol, new int [a.length]);
		
		// Prepare the decomposed matrix for substitution process
		// Given [A]{x} = {B}, let the elements of matrix {x} represents the inverse of [A].
		// Then, those of {B} are of a unit matrix
		for (int i = 0; i < a.length; ++i) {
			for (int j = 0; j < b.length; ++j)
				b[j] = (i == j ? 
						new ComplexNumber (1.0, 0.0) : new ComplexNumber (0, 0));
			
			// call substitution method
			substitute(lu, pivot, b, x[i]);
		}
		
		// return the transposed value of x
		return transpose (x);
	}
	
	/**
	 * This method returns the inverse of matrix Ai such that 
	 * 		Ai x A = In
	 * where In is the identity matrix. This method calls the standard
	 * inverse method {@link(inverse)}
	 * 
	 * @param a
	 * @return
	 * @throws MatrixException
	 */
	public static double [] [] inverse (double [] [] a) throws MatrixException
	{
		return inverse (a, new int [a.length], new double [a.length], 1.0E-18);
	}
	
	
	/**
	 * This method uses Gauss Elimination method to determine the inverse of a 
	 * square matrix [A].
	 * @param A The matrix whose inverse is desired
	 * @param pivot The transformation sequence of the rows during pivot operation
	 * @param scale The scaling factor to enhance elimination process
	 * @param tol The tolerance, which should be as low as possible to avoid division by zero
	 * @return The inverse of matrix A i.e [A]'
	 * @throws MatrixException If any error occurs, especially if a value less than tol is obtained
	 */
	public static double [] [] inverse (double [] [] a, 
			int [] pivot, double [] scale, double tol) throws MatrixException
	{
		// local variable that separates the argument from the return
		double [] [] x = new double [a.length][a.length];
		double [] b = new double [a.length];
		
		// decompose matrix a
		double [] [] lu = decompose(a, pivot, scale, tol, new int [a.length]);
		
		// Prepare the decomposed matrix for substitution process
		// Given [A]{x} = {B}, let the elements of matrix {x} represents the inverse of [A].
		// Then, those of {B} are of a unit matrix
		for (int i = 0; i < a.length; ++i) {
			for (int j = 0; j < b.length; ++j)
				b[j] = (i == j ? 1.0 : 0.0);
			
			// call substitution method
			substitute(lu, pivot, b, x[i]);
		}
		
		// return the transposed value of x
		return transpose (x);
	}
	
	/**
	 * This method uses Gauss elimination method to seperate matrix A into its corresponding
	 * lower and upper triangular matrices. The elements of the formed upper triangular matrix
	 * represent the U(i,j) from i = 0 to n and j = i to n, while the lower triangular matrix
	 * are L(i, j) from i = 1 to n for j = 0 to i.
	 * Note: The diagonal element of L are all 1's- so there's no need for representing them
	 * in a real matrix.
	 * 
	 * @param a The matrix to be decomposed into upper and lower matrices
	 * @param pivot The pivot transformation order. During decomposition, rows are supposed to
	 * be interchanged to enhance the elimination process. But instead, rows are retained but 
	 * the indices of the pivot are stored in this argument
	 * @param scale The scaling factor of the matrix. This is the maximum element in a row
	 * @param tol The tolerance, below which an element is assumed to be zero.
	 * @param multiplier Integer that tracks row interchange for computation of 
	 * the determinant of a matrix
	 * @return The decomposed matrix that has the upper and lower triangular matrix
	 * @throws MatrixException
	 */
	public static double [] [] decompose (double [] [] a, int [] pivot, double [] scale,
			double tol, int [] multiplier)	throws MatrixException
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
		for (int i = 0; i < a.length; ++i)
			for (int j = 0; j < a[i].length; ++ j)
				c[i][j] = a[i][j];
		
		// set the scale and the pivot to reflect the normal indices
		setScale (c, pivot, scale);
		
		// now begin the elimination process
		for (int k = 0; k < c.length - 1; ++k) {
			// call method setPivot
			setPivot (c, pivot, scale, k, multiplier);
			
			// check the system's tolerance before elimination
			double scaledValue = c[pivot[k]][k] / scale[pivot[k]];
			if (Math.abs(scaledValue) < tol)
				throw new MatrixException ("Scaled value " + scaledValue + 
						" is almost equal to zero: Less than the required tolerance " + tol);
			
			// operate on the lower rows
			for (int i = k + 1; i < c.length; ++i) {
				// evaluate the multiplication factor
				double factor = c [pivot[i]][k] / c [pivot[k]][k];
				
				// set the lower triangular matrix element
				c[pivot[i]][k] = factor;
				
				// eliminate pivot row from the current row
				for (int j = k + 1; j < c.length; ++j)
					c[pivot[i]][j] = c [pivot[i]][j] - factor * c [pivot[k]][j];
			}
		}
		
		return c;
	}
	
	/**
	 * This method sets the scale of the row operation during Gaussian elimination
	 * method. It seeks in each rows the column that has the highest element.
	 * @param c The matrix whose scales are to be set
	 * @param pivot	The pivot of the matrix c
	 * @param scale The values in the scale
	 */
	private static void setScale (double [] [] c, int [] pivot, double [] scale)
	{
		// iterate through the rows
		for (int i = 0; i < c.length; ++i) {
			pivot[i] = i;
			
			// set the first element of the row to scale
			scale [i] = Math.abs(c[i][0]);
			
			// choose the maximum element in the row as the scale
			for (int j = 1; j < c.length; ++j)
				if (Math.abs(c[i][j]) > scale [i])
					scale [i] = Math.abs(c[i][j]);
		}
	}
	/**
	 * This method sets the pivot for the Gauss elimination method.
	 * @param a This is the Matrix to be be decomposed by Gauss elimination method
	 * @param pivot This is the row transformation sequence that is to be set for this row operation
	 * @param scale The scaling value for the rows
	 * @param index The current Row operation that requires pivoting
	 * @param multiplier Integer that tracks row interchanger, used in computation of 
	 * the determinant
	 */
	private static void setPivot (double [] [] a, int [] pivot, double [] scale,
			int index, int [] multiplier)
	{
		// check if index is the last row of the array
		if ((index + 1) >= a.length)
			return;
		
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
		
		// check if rows have been interchanged
		if (rowPivot == index)
			return;
		
		// interchange the content of pivot at index and rowPivot
		int dummyIndex = pivot [rowPivot];
		pivot [rowPivot] = pivot[index];
		pivot [index] = dummyIndex;
		multiplier [index] = -1;;
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
	public static void substitute (ComplexNumber [] [] a, int [] pivot, 
			ComplexNumber [] b, ComplexNumber [] x)
	{
		// local variable
		int n = a.length;
		ComplexNumber sum;
		
		// begin the forward substitution
		for (int i = 0; i < n; ++i) {
			sum = b[pivot[i]];
			
			// evaluate the summation and perform the subtraction
			for (int j = 0; j < i; ++j)
				sum = sum.subtract(a[pivot[i]][j].product(b[pivot[j]]));
			
			// set the value of d for i
			b[pivot[i]] = sum;
		}
		
		// begin the backward substitution
		x[n - 1] = b[pivot[n - 1]].divide(a [pivot[n - 1]][n - 1]);
		
		// evaluate the values of matrix {x}
		for (int i = n - 2; i >= 0; --i) {
			sum = new ComplexNumber (0.0, 0.0);
			
			// compute the summation
			for (int j = i + 1; j < n; ++j)
				sum = sum.add(a [pivot[i]][j].product(x[j]));
			
			x[i] = (b[pivot[i]].subtract(sum)).divide(a[pivot[i]][i]);
		}
	}
	
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
	public static void substitute (double [] [] a, int [] pivot, double [] b, double [] x)
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
	/**
	 * This method returns the absolute value of the maximum difference 
	 * in the elements of matrix a and b. This method is useful in 
	 * determining the maximum error in matrices.
	 * @param a
	 * @param b
	 * @return The absolute value of the maximum difference in a and b
	 * @throws If the length of a is not the same as b
	 */
	public static double getMaxError (double [] a, double [] b)
	{
		if (a.length != b.length)
			throw new IllegalArgumentException ("Incompatible elements! The length " +
					"of the matrices are not the same");
		
		double error = 0;
		for (int i = 0; i < a.length; ++i)
			if (Math.abs(a[i] - b[i]) > error)
				error = Math.abs(a[i] - b[i]);
			
		return error;
	}
	
	/**
	 * This method perform multiplication of a n x n matrix by a row matrix
	 * @param a The pre-multiplier whose column number must be equal to the number 
	 * of rows of the post-multiplier. If this condition is not true, then an
	 * exception is thrown.
	 * @param b
	 * @return
	 * @throws MatrixException
	 */
	public static ComplexNumber [] product (ComplexNumber [] [] a, 
			ComplexNumber [] b)
	{
		// check validity of the arguments
		if (a == null || b == null)
			throw new MatrixException (
					"Matrix " + (a == null ? "A " : "B ") + "is empty");
		if (a[0].length != b.length)
			throw new MatrixException (
					"The columns of matrix A is not equal to rows of B");
		
		// local variables
		ComplexNumber [] c = new ComplexNumber [a.length];
		ComplexNumber sum;
		
		// perform the multiplication
		for (int i = 0; i < a.length; ++i) {
			sum  = new ComplexNumber();
			for (int j = 0; j < b.length; ++j)
				sum = sum.add(a[i][j].product(b[j]));
			c[i] = sum;
		}
		return c;
	}
	
	/**
	 * This method performs the full matrix multiplication for two matrices of 
	 * equal rows and columns
	 * @param a
	 * @param b
	 * @return
	 * @throws MatrixException If the either the number of rows of a is not equal to
	 * the number of columns of b.
	 */
	public static ComplexNumber [] [] product (ComplexNumber [] [] a, 
			ComplexNumber [] [] b)
	{
		// check validity of the arguments
		if (a == null || b == null)
			throw new MatrixException (
					"Matrix " + (a == null ? "A " : "B ") + "is empty");
		if (a[0].length != b.length)
			throw new MatrixException (
					"The columns of matrix A is not equal to rows of B");
		
		// local variables
		ComplexNumber [] [] c = new ComplexNumber [a.length] [b[0].length];
		ComplexNumber sum;
		
		for (int i = 0; i < a.length; ++i)
			for (int j = 0; j < b[i].length; ++j) {
				sum = new ComplexNumber ();
				for (int k = 0; k < a[i].length; ++k)
					sum = sum.add(a[i][k].product(b [k][j]));
				c[i][j] = sum;
			}
		return c;
	}
	
	/**
	 * This method perform multiplication of a n x n matrix by a row matrix
	 * @param a The pre-multiplier whose column number must be equal to the number 
	 * of rows of the post-multiplier. If this condition is not true, then an
	 * exception is thrown.
	 * @param b
	 * @return
	 * @throws MatrixException
	 */
	public static double [] product (double [] [] a, double [] b)
	{
		// check validity of the arguments
		if (a == null || b == null)
			throw new MatrixException (
					"Matrix " + (a == null ? "A " : "B ") + "is empty");
		if (a[0].length != b.length)
			throw new MatrixException (
					"The columns of matrix A is not equal to rows of B");
		
		// local variables
		double [] c = new double [a.length];
		
		// perform the multiplication
		for (int i = 0; i < a.length; ++i)
			for (int j = 0; j < b.length; ++j)
				c[i] = c[i] + a[i][j] * b[j];
		
		return c;
	}
	/**
	 * This method performs the full matrix multiplication for two matrices of 
	 * equal rows and columns
	 * @param a
	 * @param b
	 * @return
	 * @throws MatrixException If the either the number of rows of a is not equal to
	 * the number of columns of b.
	 */
	public static double [] [] product (double [] [] a, double [] [] b)
	{
		// check validity of the arguments
		if (a == null || b == null)
			throw new MatrixException (
					"Matrix " + (a == null ? "A " : "B ") + "is empty");
		if (a[0].length != b.length)
			throw new MatrixException (
					"The columns of matrix A is not equal to rows of B");
		
		// local variables
		double [] [] c = new double [a.length] [b[0].length];
		
		for (int i = 0; i < a.length; ++i)
			for (int j = 0; j < b[i].length; ++j)
				for (int k = 0; k < a[i].length; ++k)
					c[i][j] = c[i][j] + a[i][k] * b [k][j];
		
		return c;
	}
	
	/**
	 * This method performs the product of two 2-dimensional ExpNodes.
	 * 
	 * @param a The first ExpNode
	 * @param b The second ExpNode.
	 * 
	 * @return The product of a * b;
	 * @throws If either of the matrices are null, or the rows in matrix a 
	 * does not correspond to columns in b and vice versa
	 */
	public static ExpNode [] [] product (ExpNode [][] a, ExpNode [] [] b)
	{
		// check validity of the arguments
		if (a == null || b == null)
			throw new MatrixException (
					"Matrix " + (a == null ? "A " : "B ") + "is empty");
		if (a[0].length != b.length)
			throw new MatrixException (
					"The columns of matrix A is not equal to rows of B");
		
		// local variables
		ExpNode [] [] c = new ExpNode [a.length] [b[0].length];
		ExpNode [] p = new ExpNode [b.length];
		
		for (int i = 0; i < a.length; ++i)
			for (int j = 0; j < b[i].length; ++j) {
				for (int k = 0; k < b.length; ++k)
					p[k] = BinOpNode.make ('*', a[i][k], b [k][j]);
				
				c[i][j] = p[0];
				for (int k = 1; k < p.length; ++k)
					c[i][j] = BinOpNode.make ('+', c[i][j], p[k]);
			}
		
		return c;
	}
	
	/**
	 * This method performs the product of a 2-dimensional ExpNodes and a 
	 * column ExpNode matrices
	 * 
	 * @param a The 2-dimensional ExpNode matrix
	 * @param b The column ExpNode matrix
	 * 
	 * @return The product of a * b;
	 * @throws If either of the matrices are null, or the rows in matrix a 
	 * does not correspond to columns in b and vice versa
	 */
	public static ExpNode [] product (ExpNode [][] a, ExpNode [] b)
	{
		// check validity of the arguments
		if (a == null || b == null)
			throw new MatrixException (
					"Matrix " + (a == null ? "A " : "B ") + "is empty");
		if (a[0].length != b.length)
			throw new MatrixException (
					"The columns of matrix A is not equal to rows of B");
		
		// local variables
		ExpNode [] c = new ExpNode [a.length];
		ExpNode [] p = new ExpNode [b.length];
		
		for (int i = 0; i < a.length; ++i) {
			for (int j = 0; j < b.length; ++j)
				p[j] = BinOpNode.make ('*', a[i][j], b [j]);
			
			c[i] = p[0];
			for (int j = 1; j < p.length; ++j)
				c[i] = BinOpNode.make ('+', c[i], p[j]);
		}
		
		return c;
	}
	
	/**
	 * This method performs the product of a row matrix and a 2-dimensional
	 * ExpNode matrix
	 * 
	 * @param a The row ExpNode matrix
	 * @param b The 2-dimensional ExpNode matrix
	 * 
	 * @return The product of a * b;
	 * @throws If either of the matrices are null, or the rows in matrix a 
	 * does not correspond to columns in b and vice versa
	 */
	public static ExpNode [] product (ExpNode [] a, ExpNode [] [] b)
	{
		// check validity of the arguments
		if (a == null || b == null)
			throw new MatrixException (
					"Matrix " + (a == null ? "A " : "B ") + "is empty");
		if (a.length != b.length)
			throw new MatrixException (
					"The columns of matrix A is not equal to rows of B");
		
		// local variables
		ExpNode [] c = new ExpNode [b[0].length];
		ExpNode [] p = new ExpNode [b.length];
		
		for (int i = 0; i < a.length; ++i) {
			for (int j = 0; j < b[i].length; ++j)
				p[j] = BinOpNode.make ('*', a[j], b [j][i]);
			
			c[i] = p [0];
			for (int j = 1; j < p.length; ++j)
				c[i] = BinOpNode.make ('+', c[i], p[j]);
		}
		
		return c;
	}
	
	/**
	 * This method performs the product of a 2-dimensional int array and a 
	 * 2-dimensional ExpNode array.
	 * 
	 * @param a The argument which is the int
	 * @param b The second ExpNode.
	 * 
	 * @return The product of a * b;
	 * @throws If either of the matrices are null, or the rows in matrix a 
	 * does not correspond to columns in b and vice versa
	 */
	public static ExpNode [] [] product (int [][] a, ExpNode [] [] b)
	{
		// check validity of the arguments
		if (a == null || b == null)
			throw new MatrixException (
					"Matrix " + (a == null ? "A " : "B ") + "is empty");
		if (a[0].length != b.length)
			throw new MatrixException (
					"The columns of matrix A is not equal to rows of B");
		
		// local variables
		ExpNode [] [] c = new ExpNode [a.length] [b[0].length];
		ExpNode [] p = new ExpNode [b.length];
		
		for (int i = 0; i < a.length; ++i)
			for (int j = 0; j < b[i].length; ++j) {
				for (int k = 0; k < b.length; ++k)
					p[k] = BinOpNode.make ('*', new ConstantNode(a[i][k]), b[k][j]);
				
				c[i][j] = p[0];
				for (int k = 1; k < p.length; ++k)
					c[i][j] = BinOpNode.make ('+', c[i][j], p[k]);
			}
		
		return c;
	}
	
	/**
	 * This method performs the product of a row int array and a 
	 * 2-dimensional ExpNode array. The returned matrix is a 
	 * row ExpNode array
	 * 
	 * @param a The argument which is the int
	 * @param b The second ExpNode.
	 * 
	 * @return The product of a * b;
	 * @throws If either of the matrices are null, or the rows in matrix a 
	 * does not correspond to columns in b and vice versa
	 */
	public static ExpNode [] product (int [] a, ExpNode [] [] b)
	{
		// check validity of the arguments
		if (a == null || b == null)
			throw new MatrixException (
					"Matrix " + (a == null ? "A " : "B ") + "is empty");
		if (a.length != b.length)
			throw new MatrixException (
					"The columns of matrix A is not equal to rows of B");
		
		// local variables
		ExpNode [] c = new ExpNode [b[0].length];
		ExpNode [] p = new ExpNode [b.length];
		
		for (int i = 0; i < b[0].length; ++i) {
			for (int j = 0; j < b.length; ++j)
				p[j] = BinOpNode.make ('*', new ConstantNode(a[j]), b [j][i]);
			
			c[i] = p[0];
			for (int j = 1; j < p.length; ++j)
				c[i] = BinOpNode.make ('+', c[i], p[j]);
		}
		
		return c;
	}
	
	/**
	 * This method returns the product of a 2-dimensional ExpNode array and 
	 * a 2-dimensional int array.
	 * 
	 * @param a The first ExpNode
	 * @param b The second argument which is an int array.
	 * 
	 * @return The product of a * b;
	 * @throws If either of the matrices are null, or the rows in matrix a 
	 * does not correspond to columns in b and vice versa
	 */
	public static ExpNode [] [] product (ExpNode [][] a, int [] [] b)
	{
		// check validity of the arguments
		if (a == null || b == null)
			throw new MatrixException (
					"Matrix " + (a == null ? "A " : "B ") + "is empty");
		if (a[0].length != b.length)
			throw new MatrixException (
					"The columns of matrix A is not equal to rows of B");
		
		// local variables
		ExpNode [] [] c = new ExpNode [a.length] [b[0].length];
		ExpNode [] p = new ExpNode [b.length];
		
		for (int i = 0; i < c.length; ++i)
			for (int j = 0; j < b[i].length; ++j) {
				for (int k = 0; k < b.length; ++k) 
					p[k] = BinOpNode.make ('*', a[i][k], new ConstantNode(b[k][j]));
				
				c[i][j] = p[0];
				for (int k = 1; k < p.length; ++k)
					c[i][j] = BinOpNode.make ('+', c[i][j], p[k]);
			}
		
		return c;
	}
	
	/**
	 * This method performs the product of a 2-dimensional int array and 
	 * a column ExpNode. The returned matrix is a column matrix
	 * 
	 * @param a The int array ExpNode
	 * @param b The column ExpNode matrix.
	 * 
	 * @return The product of a * b;
	 * @throws If either of the matrices are null, or the columns in matrix a 
	 * does not correspond to rows in b.
	 */
	public static ExpNode [] product (int [][] a, ExpNode [] b)
	{
		// check validity of the arguments
		if (a == null || b == null)
			throw new MatrixException (
					"Matrix " + (a == null ? "A " : "B ") + "is empty");
		if (a[0].length != b.length)
			throw new MatrixException (
					"The columns of matrix A is not equal to rows of B");
		
		// local variables
		ExpNode [] c = new ExpNode [a.length];
		ExpNode [] p = new ExpNode [b.length];
		
		for (int i = 0; i < a.length; ++i) {
			for (int j = 0; j < b.length; ++j)
				p[j] = BinOpNode.make ('*', new ConstantNode(a[i][j]), b [j]);
			
			c[i] = p[0];
			for (int j = 1; j < p.length; ++j)
				c[i] = BinOpNode.make ('+', c[i], p[j]);
		}
		
		return c;
	}
	
	/**
	 * This method performs the product of a row ExpNode matrix and 
	 * a 2-dimensional int array. The returned ExpNode is a row matrix.
	 * 
	 * @param a The row Matrix ExpNode
	 * @param b The 2-dimensional int array
	 * 
	 * @return The product of a * b;
	 * @throws If the columns of a does not correspond to the rows of b
	 */
	public static ExpNode [] product (ExpNode [] a, int [] [] b)
	{
		// check validity of the arguments
		if (a == null || b == null)
			throw new MatrixException (
					"Matrix " + (a == null ? "A " : "B ") + "is empty");
		if (a.length != b.length)
			throw new MatrixException (
					"The columns of matrix A is not equal to rows of B");
		
		// local variables
		ExpNode [] c = new ExpNode [b[0].length];
		ExpNode [] p = new ExpNode [b.length];
		
		for (int i = 0; i < b[0].length; ++i) {
			for (int j = 0; j < b.length; ++j)
				p[j] = BinOpNode.make ('*', a[j], new ConstantNode (b [j][i]));
			
			c[i] = p[0];
			for (int j = 1; j < p.length; ++j)
				c[i] = BinOpNode.make ('+', c[i], p[j]);
		}
		
		return c;
	}
	
	/**
	 * This method returns an ExpNode that is a product of a row ExpNode matrix
	 * and a column ExpNode matrix.
	 * 
	 * @param a The first ExpNode
	 * @param b The second ExpNode.
	 * 
	 * @return The product of a * b;
	 * @throws If either of the matrices are null, length of matrix a and b 
	 * are not equal
	 */
	public static ExpNode product (ExpNode [] a, ExpNode [] b)
	{
		// check validity of the arguments
		if (a == null || b == null)
			throw new MatrixException (
					"Matrix " + (a == null ? "A " : "B ") + "is empty");
		if (a.length != b.length)
			throw new MatrixException (
					"The columns of matrix A is not equal to rows of B");
		
		// local variables
		ExpNode c = BinOpNode.make ('*', a[0], b[0]);
		
		for (int i = 1; i < a.length; ++i)
			c = BinOpNode.make ('+', c, BinOpNode.make ('*', a[i], b[i]));
		
		return c;
	}
	
	/**
	 * This method returns an ExpNode that is a product of a row ExpNode matrix
	 * and a column int matrix.
	 * 
	 * @param a The ExpNode array
	 * @param b The int array.
	 * 
	 * @return The product of a * b;
	 * @throws If either of the matrices are null, length of matrix a and b 
	 * are not equal
	 */
	public static ExpNode product (ExpNode [] a, int [] b)
	{
		// check validity of the arguments
		if (a == null || b == null)
			throw new MatrixException (
					"Matrix " + (a == null ? "A " : "B ") + "is empty");
		if (a.length != b.length)
			throw new MatrixException (
					"The columns of matrix A is not equal to rows of B");
		
		// local variables
		ExpNode c = BinOpNode.make ('*', a[0], new ConstantNode (b[0]));
		
		for (int i = 1; i < a.length; ++i)
			c = BinOpNode.make ('+', c, BinOpNode.make ('*', 
					a[i], new ConstantNode (b[i])));
		
		return c;
	}
	
	/**
	 * This method returns an ExpNode that is a product of a row int matrix
	 * and a column ExpNode matrix.
	 * 
	 * @param a The int array
	 * @param b The ExpNode array.
	 * 
	 * @return The product of a * b;
	 * @throws If either of the matrices are null, length of matrix a and b 
	 * are not equal
	 */
	public static ExpNode product (int [] a, ExpNode [] b)
	{
		// check validity of the arguments
		if (a == null || b == null)
			throw new MatrixException (
					"Matrix " + (a == null ? "A " : "B ") + "is empty");
		if (a.length != b.length)
			throw new MatrixException (
					"The columns of matrix A is not equal to rows of B");
		
		// local variables
		ExpNode c = BinOpNode.make ('*', new ConstantNode (a[0]), b[0]);
		
		for (int i = 1; i < a.length; ++i)
			c = BinOpNode.make ('+', c, BinOpNode.make ('*', 
					new ConstantNode (a[i]), b[i]));
		
		return c;
	}
	
	/**
	 * This method performs matrix subtraction by subtracting the values of b from a i.e.
	 * subtract (a, b) = c
	 * where
	 *  [C] = [A] - [B]
	 * @param a
	 * @param c
	 * @return
	 */
	public static double [] [] subtract (double [] [] a, double [] [] b)
	{
		// local variable
		double [] [] c = new double [a.length] [b[0].length];
		
		for (int i = 0; i < a.length; ++i)
			for (int j = 0; j < a[0].length; ++j)
				c[i][j] = a[i][j] - b[i][j];
		
		return c;
	}

	/**
	 * This method performs matrix subtraction by subtracting the values of b from a i.e.
	 * subtract (a, b) = c
	 * where
	 *  {C} = {A} - {B}
	 * @param a
	 * @param c
	 * @return
	 */
	public static double [] subtract (double [] a, double [] b)
	{
		// local variable
		double [] c = new double [a.length];
		
		for (int i = 0; i < a.length; ++i)
			c[i] = a[i] - b[i];
		
		return c;
	}
	
	/**
	 * This method computes the scalar product of a matrix
	 * @param k
	 * @param a
	 * @return
	 */
	public static double [] [] product (double k, double [] [] a)
	{
		double [] [] c = new double [a.length][a[0].length];
		for (int i = 0; i < c.length; ++i)
			for (int j = 0; j < a[i].length; ++j)
				c[i][j] = k * a[i][j];
		return c;
	}
	
	/**
	 * This method transposes the elements of the square matrix a.
	 * @param a
	 * @return
	 */
	public static ComplexNumber [] [] transpose (ComplexNumber [] [] a)
	{
		ComplexNumber [] [] t = new ComplexNumber [a[0].length] [a.length];
		
		for (int i = 0; i < a.length; ++i)
			for (int j = 0; j < a[i].length; ++j)
				t [i] [j] = a [j][i];
		
		return t;
	}
	
	/**
	 * This method transposes a row matrix to a column matrix
	 * @param a
	 * @return
	 */
	public static ComplexNumber [] [] transpose (ComplexNumber [] a)
	{
		ComplexNumber [] [] t = new ComplexNumber [a.length][1];
		
		for (int i = 0; i < a.length; ++i)
			t [i][0] = a[i];
		
		return t;
	}
	
	/**
	 * This method transposes the elements of the square matrix a.
	 * @param a
	 * @return
	 */
	public static double [] [] transpose (double [] [] a)
	{
		double [] [] t = new double [a[0].length] [a.length];
		
		for (int i = 0; i < a.length; ++i)
			for (int j = 0; j < a[i].length; ++j)
				t [i] [j] = a [j][i];
		
		return t;
	}
	
	/**
	 * This method transposes a row matrix to a column matrix
	 * @param a
	 * @return
	 */
	public static double [] [] transpose (double [] a)
	{
		double [] [] t = new double [a.length][1];
		
		for (int i = 0; i < a.length; ++i)
			t [i][0] = a[i];
		
		return t;
	}
	
	
	/**
	 * This method computes the Lyapunov exponents for a continuous nonlinear
	 * system whose state functions are described by ExpNode fx.
	 * 
	 * @param f The state function of the continuous nonlinear system
	 * @param x The state variables of the system
	 * @param iterate The number of iterates that should be performed
	 * @param variables The variables that are defined in fx.
	 * @param values The initial values of the state variables above
	 * @param dt The step value of the rungeKutta iteration
	 * 
	 * @return The Lyapunove exponents.
	 */
	public static double [] getLyapunoveExpCont (ExpNode [] f, String [] x, 
			int iterate, String [] variables, double [] values, double dt)
	{
		// local variables
		int n = x.length;
		double [] h = new double [n];
		String [] var = Matrix.getVariationVar(variables, x);
		double [] val = new double [values.length];
		double [] val1 = new double [0];
		
		// copy content of values into val
		for (int i = 0; i < values.length; ++i) val [i] = values[i];
		
		// for the array of the independent coordinates x
		int [] xIndex = ExpNode.getXIndex(
				getVariationIndependentVar(x), var);
		
		// form the ExpNode for the first variational equation
		ExpNode [] fx = Matrix.getVariationEq(f, x);
		
		// perform a hundred iteration
		for (int i = 0; i < 1000; ++i)
			val = Algebra.rungeKutta(f, xIndex, var, val, dt);
		
		double [] [] w = getIdentity (x.length);
		double [] [] z, Df;
		
		// perform Lyapunov computation via Gram-Schmidt Orthognalization
		for (int i = 0; i < iterate; ++i) {
			// reset the variational variables
			val1 = Matrix.getFirstVariationInitVal(val, x);
			
			// get the time-1 iterate
			for (int t = 0; t < ((int) 1.0 / dt); ++t)
				val1 = Algebra.rungeKutta(fx, xIndex, var, val1, dt);
			
			// copy the iterated values of val1 into val
			for (int k = 0; k < n; ++k)
				val[xIndex[k]] = val1[xIndex[k]];
			
			// get the ellipsoid
			Df = new double [n] [n];
			for (int index = 0; index < n * n; ++index)
				Df [index / n] [index % n] = val1 [xIndex[n + index]];
			
			z = Matrix.product(Df, w);
			
			// orthogonalize z
			w = getOrthogonalVector (z);
			
			// get the expansion and convert to lyapunov exponents
			double [] len = getModulus (w);
			for (int j = 0; j < len.length; ++j)
				h [j] += Math.log(len[j]);
			
			// normalize the vector
			normalize (w);
		}
		
		// divide the exponents by number of iterations
		for (int i = 0; i < n; ++i)
			h[i] /= iterate;
		return h;
	}
	
	
	/**
	 * This method computes the Lyapunov exponents for a discrete nonlinear
	 * system whose state functions are described by ExpNode fx.
	 * 
	 * @param fx The state function of the continuous nonlinear system
	 * @param x The state variables of the discrete system
	 * @param iterate The number of iterates that should be performed
	 * @param var The variables that are defined in fx.
	 * @param values The initial values of the state variables above
	 * @return The Lyapunove exponents.
	 */
	public static double [] getLyapunoveExpDiscrete (ExpNode [] fx, String [] x, 
			int iterate, String [] var, double [] values)
	{
		// local variables
		int n = x.length;
		double [] h = new double [n];
		double [] val = new double [values.length];
		double [] temp = new double [n];
		int [] xIndex = ExpNode.getXIndex(x, var);
		ExpNode [] [] Df = Matrix.jacobian(fx, x);
		
		// copy content of values into val
		for (int i = 0; i < val.length; ++i)
			val [i] = values [i];
		
		// perform a hundred iteration
		for (int i = 0; i < 100; ++i) {
			temp = Matrix.values(fx, var, val);
			for (int j = 0; j < n; ++j)
				val[xIndex[j]] = temp[j];
		}
		
		double [] [] w = getIdentity (n);
		double [] [] z, df;
		
		// perform the number of iterations
		for (int i = 0; i < iterate; ++i) {
			// get the new iterate
			temp = Matrix.values(fx, var, val);
			for (int j = 0; j < n; ++j)
				val[xIndex[j]] = temp[j];
			
			// get the ellipsoid
			df = Matrix.values(Df, var, val);
			z = product (df, w);
			
			// orthogonalize z
			w = getOrthogonalVector (z);
			
			// get the expansion and convert to lyapunov exponents
			double [] len = getModulus (w);
			for (int j = 0; j < len.length; ++j)
				h [j] += Math.log(len[j]);
			
			// normalize the vector
			normalize (w);
		}
		
		// divide the exponents by number of iterations
		for (int i = 0; i < n; ++i)
			h[i] /= iterate;
		return h;
	}
	
	/**
	 * This method return the orthogonal matrix of the vector basis z using
	 * the Gram-Schmidth orthogonalization scheme.
	 * @param z
	 * @return
	 */
	public static double [] [] getOrthogonalVector (double [] [] z)
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
	public static void normalize (double [] [] y)
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
	public static double [] getModulus (double [] [] y)
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
	 * This method returns the identity matrix for n X n matrix
	 * @param n
	 * @return
	 */
	public static double [] [] getIdentity (int n)
	{
		// local variable
		double [] [] in = new double [n] [n];
		
		for (int i = 0; i < n; ++i)
			in [i] [i] = 1.0;
		
		return in;
	}
	
	/**
	 * This method returns the derivative of the special section Pi that is
	 * used in mapping Poicare section
	 * @param x The independent coordinates
	 * @param sectionIndex The Poicare section
	 * @return derivative of the embedded section mapping
	 */
	public static int [] [] Dhx (String [] x, int sectionIndex)
	{
		int n = x.length;
		int [] [] dh = new int [n - 1][n];
		
		for (int i = 0, j = 0; i < n; ++i)
			if (i != sectionIndex)
				dh[j++][i] = 1;
		
		return dh;
	}
	
	/**
	 * This method returns the derivative of the special section Pi that is
	 * used in mapping Poicare section
	 * 
	 * @param x The independent coordinates
	 * @param sectionIndex The Poicare section
	 * @return derivative of the embedded section mapping
	 */
	public static int [] Dqx (String [] x, int sectionIndex)
	{
		int n = x.length;
		int [] dq = new int [n];
		
		// set the index of section to 1.
		dq [sectionIndex] = 1;
		
		return dq;
	}
	
	/**
	 * This method returns the derivative of the return map of the embedded 
	 * section used in Poincares section
	 * @param x The independent coordinates
	 * @param sectionIndex The Poicare section
	 * @return derivative of the embedded section mapping
	 */
	public static int [] [] Dhu (String [] x, int sectionIndex)
	{
		int n = x.length;
		int [] [] dhu = new int [n][n - 1];
		
		for (int i = 0, j = 0; i < n; ++i)
			if (i != sectionIndex)
				dhu[i][j++] = 1;
		
		return dhu;
	}
	
	/**
	 * This method returns the embedded matrix used in the Poincare section 
	 * mapping
	 * @param x
	 * @param sectionIndex
	 * @return
	 */
	public static String [] hx (String [] x, int sectionIndex)
	{
		int n = x.length;
		String [] h = new String [n - 1];
		for (int i = 0, j = 0; i < n; ++i)
			if (i != sectionIndex)
				h[j++] = x[i];
		return h;
	}
	
	/**
	 * This method allows displaying of a one-dimensional matrix on the console
	 * @param a
	 */
	public static void display (ComplexNumber [] a)
	{
		for (int i = 0; i < a.length; ++i)
			System.out.print(a[i] + "   ");
		System.out.println();	
	}
	
	/**
	 * This method displays a two-dimensional matrix on the console
	 * @param a
	 */
	public static void display (ComplexNumber [] [] a)
	{
		for (int i = 0; i < a.length; ++i) {
			for (int j = 0; j < a[i].length; ++j)
				System.out.print(a[i][j] + "   ");
			System.out.println();
		}
		System.out.println();
	}
	
	/**
	 * This method allows displaying of a one-dimensional matrix on the console
	 * @param a
	 */
	public static void display (double [] a)
	{
		for (int i = 0; i < a.length; ++i)
			System.out.print(a[i] + "   ");
		System.out.println();
	}
	
	/**
	 * This method displays a two-dimensional matrix on the console
	 * @param a
	 */
	public static void display (double [] [] a)
	{
		for (int i = 0; i < a.length; ++i) {
			for (int j = 0; j < a[i].length; ++j)
				System.out.print(a[i][j] + "   ");
			System.out.println();
		}
		System.out.println();
	}
	
	/**
	 * This method displays a String array
	 * @param v
	 */
	public static void display (String [] v) {
		for (int i = 0; i < v.length; ++i)
			System.out.print(v[i] + "   ");
		System.out.println();
	}
	
	/**
	 * This method displays a 2-dimensional String array
	 * @param v
	 */
	public static void display (String [] [] v) {
		for (int i = 0; i < v.length; ++i)
			display (v[i]);
	}
	
	/**
	 * This method displays an ExpNode
	 * @param a
	 */
	public static void displayInfix (ExpNode [] [] a)
	{
		for (int i = 0; i < a[0].length; ++i) {
			System.out.println("Column " + (i + 1));
			for (int j = 0; j < a.length; ++j)
				System.out.println("\t" + a[j][i].printInfix());
		}
	}
	
	/**
	 * This method displays an ExpNode as a column matrix
	 * @param a
	 */
	public static void display (ExpNode [] a)
	{
		for (ExpNode f : a)
			System.out.println ("\t" + f.printInfix());
	}
}
