package bifurcation;

import java.awt.Point;

import math.tools.ComplexNumber;
import math.tools.Matrix;
import math.tools.MatrixException;


public class ContinunityTest 
{
	// instance variables
	private double g;
	private final int TANGENT = 1;
	private final int PERIOD = 2;
	private final int NEIMARK = 3;
	private final int HOPF = 4;
	
	private class BifurcType {
		double y,	// the fixed point
			   t,		// period of the limit cycle
			   k,		// the k-value at the bifurcation
			   type,	// the bifurcation type
			   h;		// the bifurcation coefficient
		
		public BifurcType (int type, double y, double t, double k, double h)
		{
			this.type = type;
			this.y = y;
			this.t = t;
			this.k = k;
			this.h = h;
		}
		
		public void print () {
			System.out.println ("Type = " + type + ", Fixed Point = " + y +
					", Period = " + t +	", k-value = " + k + ", h = " + h);
		}
	}
	
	// constructor for limit cycle tes
	public ContinunityTest () {
		java.util.ArrayList<BifurcType> bifurcType = new java.util.ArrayList<>();
		
		// prompt the user for inputs
		java.util.Scanner scanner = new java.util.Scanner(System.in);
		scanner.useDelimiter("\r\n");
		
		try {
			// prompt for entry of the g-value at which the bifurcation points
			// will be evaluated;
			System.out.println("Enter the g-value: ");
			g = 1.6; //scanner.nextDouble();
			
			// start the k-values from the saddle node bifurcation
			// till k = 1
			double dt = 1.0E-2;
			double k = 0.0001 + 1.0 / g;
			double error = 1.0;
			while (k < 0.845) {
				// intialize a new start point on the poincare section
				double x = 4.102;
				double [] x1 = {x, 0.0, 1.0, 0.0, 0.0, 1.0};
				x1 = rungeKutta (x1, g, k, dt);	
				double [] x2 = rungeKutta (x1, g, k, dt);
				
				// get the next poincare crossing
				int count = 1;
				while (!(x1[1] <= 0.0 && x2[1] > 0.0)) {
					// copy content of x2 into x1
					x1 = x2;
					x2 = rungeKutta (x1, g, k, dt);
					++count;
				}
				
				// iterative Newton method
				error = Math.abs(x1[0] - x);
				double [] xn = {x, dt * count};
				double [] xn1;
				double [] f = {x1[0] - x, 0.0};
				double [] [] Df = {{x1[2] - 1.0, Math.tanh(g * x)}, {x1[4], x}};
				
				while (error > 1.0E-6) {
					// evaluate current xn, f and Df
					xn1 = Matrix.subtract(xn, Matrix.product(Matrix.inverse(Df), f));
					
					// iterate starting the newly computed value of x
					x1 = new double [] {xn1[0], 0.0, 1.0, 0.0, 0.0, 1.0};
					x1 = rungeKutta (x1, g, k, dt);	
					x2 = rungeKutta (x1, g, k, dt);
					
					// get the next poincare crossing
					while (!(x1[1] <= 0.0 && x2[1] > 0.0) 
							&& Math.abs(x1[0] - x2[0]) > 1.0E-6) {
						// copy content of x2 into x1
						x1 = x2;
						x2 = rungeKutta (x1, g, k, dt);
					}
					
					// confirm that the orbit has not gone to the equilibrium point
					if (Math.abs(x1[0] - x2[0]) < 1.0E-6)
						break;
					
					// evaluate the current error in xn
					//error = Matrix.getMaxError(xn1, xn);
					error = Math.abs(xn1[0] - xn[0]);
					xn = xn1;//new double [] {xn1[0], dt * count};
					
					// prepare variables for the next iterate
					f = new double [] {x1[0] - xn[0], x1[1]};
					Df = new double [] [] {
							{x1[2] - 1.0, -x1[1] + Math.tanh(g * x1[0])},
							{x1[4], x1[0] - k * x1[1]}};
				}
				
				// check the eigenvalues of x1 to see if there is a bifurcation
				ComplexNumber [] e = Matrix.getEigenValues(new double [][] {
						{x1[2], x1[3]}, {x1[4], x1[5]}});
				error = 1.0E-5;
				for (ComplexNumber h : e)
					if (Math.abs(h.getRealNumber() - 1.0) <= error && 
							Math.abs(h.getImaginaryNumber()) <= error)
						bifurcType.add(new BifurcType (TANGENT, xn[0], xn[1], k, 
								h.getRealNumber()));
					else if (Math.abs(h.getRealNumber() + 1.0) <= error &&
							Math.abs(h.getImaginaryNumber()) <= error)
						bifurcType.add(new BifurcType (PERIOD, xn[0], xn[1], k,
								h.getRealNumber()));
					else if (Math.abs(h.getRealNumber()) <= error &&
							Math.abs(h.getImaginaryNumber()) > error)
						bifurcType.add(new BifurcType (HOPF, xn[0], xn[1], k,
								h.getImaginaryNumber()));
					else if (Math.abs(h.getModulus() - 1.0) <= error)
						bifurcType.add(new BifurcType (NEIMARK, xn[0], xn[1], k,
								h.getAngleTheta ()));
				
				k += 0.0001;
				System.out.println("k = " + k);
			}
			
			// begin the continuity test
			for (BifurcType b : bifurcType) {
				// local variables
				double [] x1, x2, f, xn1;
				double [] [] Df;
				double [] xn = new double [] {b.y, b.t, b.k};
				
				System.out.println ("\nAnalysis at \n" +
						"k = " + b.k + ", u = " + b.h + "\ng\tk");
				
				// iterate for continuity
				for (double gg = g + 0.0001; gg <= 2.0; gg += 0.0001) {
					error = 1.0;
					
					// evaluate the new solution
					while (error > 1.0E-6) {
						// at the new value, re-initialize x1 and
						// iterate until next poincare crossing
						x1 = rungeKuttafc (new double [] {
								xn[0], 0.0, 1.0, 0.0, 0.0, 1.0, 
								0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
								0.0, 0.0, 0.0, 0.0, 0.0}, gg, xn[2], dt);
						x2 = rungeKuttafc (x1, gg, xn[2], dt);
						while (!(x1[1] <= 0.0 && x2[1] > 0.0)
								&& Math.abs(x1[0] - x2[0]) > 1.0E-6) {
							x1 = x2;
							x2 = rungeKuttafc (x1, gg, xn[2], dt);
						}
						
						// evaluate the new coordinate
						f = new double [] {x1[0] - xn[0], x1[1], getXu(x1, 1.0)};
						Df = getDffc (x1, gg, xn[2]);
						xn1 = Matrix.subtract(xn, Matrix.product(Matrix.inverse(Df), f));
						
						// evaluate the new error
						error = Math.abs(xn1[2] - xn[2]);
						xn = xn1;
					}
					
					// print the value of xn
					System.out.printf ("%f4\t%f4\n", gg, xn[2]);
				}
			}
			
		} catch (MatrixException e) {
			e.printStackTrace();
		} catch (NumberFormatException e) {
			e.printStackTrace();
		}
	}
	
	// continuity determinant
	private double getXu (double [] x, double u) throws MatrixException
	{
		return Matrix.getDeterminant(new double [] [] {
				{x[2] - u, x[3]},
				{x[4], x[5] - u}
		});
	}
	
	// continuity differential
	private double [] [] getDffc (double [] x, double g, double k)
	{
		// local variable
		double f1 = -x[1] + Math.tanh(g * x[0]);
		double f2 = x[0] - k * x [1];
		
		return new double [] [] {
				{x[2] - 1, f1, 0.0},
				{x[4], f2, -x[1]},
				{-x[8] - x[11], -f1 * (x[8] + x[11]) - f2 * (x[10] + x[13]),
					-x[6] * (x[8] + x[11]) - x[7] * (x[10] + x[9])}
		};
	}
	
	// the continuity function for limit cycle
	private double [] ffc (double [] v, double g, double k)
	{
		// local variables
 		double [] vv = new double [v.length];
 		double sech = 1.0 / Math.cosh(g * v[0]);
 		
 		vv[0] = -v[1] + Math.tanh(g * v[0]);
 		vv[1] = v[0] - k * v[1];
 		
 		
 		for (int i = 0; i < 2; ++i) {
 			vv[i + 2] = v[i + 2] * g * sq(sech) - v[i + 4];
 			vv[i + 4] = v[i + 2] - k * v[i + 4];
 		}
 		
 		vv[6] = g * sq(sech) * v[6] - v[7];
 		vv[7] = v[6] - k * v[7] - v[1];
 		
 		vv[8] = g * sq(sech) * v[8] - v[9] 
 				- 2 * sq(g) * Math.tanh(g * v[0]) * sq (sech) * sq(v[2]);
 		vv[9] = v[8] - k * v[9];
 		
 		vv[10] = g * sq(sech) * v[10] - v[11] 
 				- 2 * sq(g) * Math.tanh(g * v[0]) * sq (sech) * v[2] * v[3];
 		vv[11] = v[10] - k * v[11];
 		
 		vv[12] = g * sq(sech) * v[12] - v[13] 
 				- 2 * sq(g) * Math.tanh(g * v[0]) * sq(sech) * sq(v[3]);
 		vv[13] = v[12] - k * v[13];
 		
 		vv[14] = g * sq(sech) * v[14] - v[15] 
 				- 2 * sq(g) * Math.tanh(g * v[0]) * sq(sech) * v[2] * v[7];
 		vv[15] = v[14] - k * v[15] - v[4];
 		
 		vv[16] = g * sq(sech) * v[16] - v[17]
 				- 2 * sq(g) * Math.tanh(g * v[0]) * v[3] * v[6];
 		vv[17] = v[16] - k * v[17] - v[5];
 		
 		return vv;
	}
	
	// utility method sq
	private double sq (double x) {return x * x; }
	
	// Runge Kutta fourth order numerical integeration
 	private double [] rungeKuttafc (double [] v, double g, double k, double dt)
 	{
 		int l = v.length;
 		double [] c1, c2, c3, c4;
 		
 		// initialize the intermediate steps
 		c1 = new double [l];
 		c2 = new double [l];
 		c3 = new double [l];
 		c4 = new double [l];
 		
 		c1 = ffc(v, g, k);
 		
 		for (int i = 0; i < l; ++i)
 			c2[i] = v[i] + dt * c1[i] / 2;
 		c2 = ffc(c2, g, k);
 		
 		for (int i = 0; i < l; ++i)
 			c3[i] = v[i] + dt * c2[i] / 2;
 		c3 = ffc(c3, g, k);
 		
 		for (int i = 0; i < l; ++i)
 			c4[i] = v[i] + dt * c3[i];
 		c4 = ffc(c4, g, k);
 		
 		for (int i = 0; i < l; ++i)
 			c1[i] = v[i] + dt * (c1[i] + 2 * (c2[i] + c3[i]) + c4[i]) / 6.0;
 		
 		return c1;
 	}

 	// Runge Kutta fourth order numerical integeration
 	private double [] rungeKutta (double [] v, double g, double k, double dt)
 	{
 		int l = v.length;
 		double [] c1, c2, c3, c4;
 		
 		// initialize the intermediate steps
 		c1 = new double [l];
 		c2 = new double [l];
 		c3 = new double [l];
 		c4 = new double [l];
 		
 		c1 = f(v, g, k);
 		
 		for (int i = 0; i < l; ++i)
 			c2[i] = v[i] + dt * c1[i] / 2;
 		c2 = f(c2, g, k);
 		
 		for (int i = 0; i < l; ++i)
 			c3[i] = v[i] + dt * c2[i] / 2;
 		c3 = f(c3, g, k);
 		
 		for (int i = 0; i < l; ++i)
 			c4[i] = v[i] + dt * c3[i];
 		c4 = f(c4, g, k);
 		
 		for (int i = 0; i < l; ++i)
 			c1[i] = v[i] + dt * (c1[i] + 2 * (c2[i] + c3[i]) + c4[i]) / 6.0;
 		
 		return c1;
 	}
 	
 	// method f(x)
 	private double [] f (double [] v, double g, double k)
 	{
 		// local variables
 		double [] vv = new double [v.length];
 		double sech = 1.0 / Math.cosh (g * v[0]);
 		
 		vv[0] = -v[1] + Math.tanh(g * v[0]);
 		vv[1] = v[0] - k * v[1];
 		
 		
 		for (int i = 0; i < 2; ++i) {
 			vv[i + 2] = v[i + 2] * g * sech * sech - v[i + 4];
 			vv[i + 4] = v[i + 2] - k * v[i + 4];
 		}
 		
 		return vv;
 	}
	
	/*
	// constructor
	public ContinunityTest () {
		double [] v0 = new double [2];
		
		// prompt the user for inputs
		java.util.Scanner scanner = new java.util.Scanner(System.in);
		scanner.useDelimiter("\r\n");
		
		// repeat process
		while (true) {
			try {
				System.out.println("Enter the parameters as g, k: ");
				//String [] param = scanner.next().split(",");
				g = 0.1; //Double.parseDouble(param[0]);
				double k = 0.1; //Double.parseDouble(param[1]);
				
				System.out.println("Enter the initial value of (x, y): ");
				// String [] initial = scanner.next().split(",");
				v0[0] = 1.8; //Double.parseDouble(initial[0]);
				v0[1] = -0.25; //Double.parseDouble(initial[1]);
				
				// call newton's shooting method to find the stable fixed point
				v0 = shootForFixedPoint (v0, k, 1.0E-6);
				
				// evaluate the jacobian at the fixed point
				double [] [] df = Df(v0, k);
				ComplexNumber [] lambda = Matrix.getEigenValues(df);
				
				// display the fixed point and its corresponding eigen value
				display(v0);
				for (ComplexNumber c : lambda)
					System.out.println(c);
				
				// print the continunity values of g and k
				System.out.println("|  g    |    k     |");
				// form a 3-element v0
				//0.5;k=2.0;
				double h = lambda[0].getImaginaryNumber();
				
				double [] vv = new double [] {v0[0], v0[1], k, h};
				for ( ; g <= 2.0; g += 0.01) {
					vv = hBifurcationContinunity (vv, 1.0E-6);
					System.out.printf ("%f4\t%f4\t%f4\n", g, vv[2], vv[3]);
				}
				System.exit(0);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}	// end constructor
	*/
	
	// shooting continunity method for hopf bifurcation
	private double [] hBifurcationContinunity (double [] v, double error)
	{
		// local variable
		double [] v0 = new double [v.length];
		for (int i = 0; i < v.length; ++i)
			v0[i] = v[i];
		
		double [] fx = new double [v0.length];
		double [][] Dfx;
		int counter = 0;
		
		try {
			// perform the iteration
			while (++counter < 100) {
				// evaluate f(x), Df(x) and new v
				fx = fh(v0);
				Dfx = Dfh(v0);
				fx = Matrix.subtract(v0, Matrix.product(Matrix.inverse(Dfx), fx));
				
				// evaluate the percentage error in v1
				if (Matrix.getMaxError(fx, v0) < error)
					return fx;
				else
					v0 = fx;
			}
			throw new MatrixException ("Failure to converge");
		} catch (MatrixException e) {
			// TODO Auto-generated catch block1
			System.out.println(e.getMessage());
		}
		return v0;
	}
	
	// function for the hopf continunity
	private double [] fh (double [] x)
	{
		double sech = 1.0 / Math.cosh(g * x[0]);
		
		return new double [] {
				-x[1] + Math.tanh(g * x[0]),
				x[0] - x[2] * x[1],
				-x[3] * x[3] - g * x[2] * sech * sech + 1,
				x[3] * (x[2] - g * sech * sech)
		};
	}
	
	// DF for hopf bifurcation continunity
	private double [] [] Dfh (double [] x)
	{
		double sech = 1.0 / Math.cosh(g * x[0]);
		
		return new double [] [] {
				{g * sech * sech, -1.0, 0.0, 0.0},
				{1.0, -x[2], -x[1], 0.0},
				{2 * g * g * x[2] * Math.tanh(g * x[0]) * sech * sech, 0, 
					g * sech * sech, -2.0 * x[3]},
				{2 * x[3] * g * Math.tanh(g * x[0]) * sech * sech, 0, x[3],
						x[2] - g * sech * sech
				}
				
		};
	}
	
	// continunity shooting method
	private double [] snBifurcationContinunity (double [] v, double u, double error)
	{
		// local variable
		double [] v0 = new double [v.length];
		for (int i = 0; i < v.length; ++i)
			v0[i] = v[i];
		
		double [] fx = new double [v0.length];
		double [][] Dfx;
		int counter = 0;
		
		try {
			// perform the iteration
			while (++counter < 100) {
				// evaluate f(x), Df(x) and new v
				fx = fc(v0, u);
				Dfx = Dfc(v0, u);
				fx = Matrix.subtract(v0, Matrix.product(Matrix.inverse(Dfx), fx));
				
				// evaluate the percentage error in v1
				if (Matrix.getMaxError(fx, v0) < error)
					return fx;
				else
					v0 = fx;
			}
			throw new MatrixException ("Failure to converge");
		} catch (MatrixException e) {
			// TODO Auto-generated catch block
			System.out.println(e.getMessage());
		}
		return v0;
	}
	
	// method for evaluating f(x, k) for bifurcation continunity
	private double [] fc (double [] x, double u) throws MatrixException
	{
		return new double [] {
				-x[1] + Math.tanh(g * x[0]),
				x[0] - x[2] * x[1],
				1.0 - (x[2] + u) * (g * Math.pow(Math.cosh(g * x[0]), -2) - u)
				/*Matrix.getDeterminant(new double [] [] {
						{g * Math.pow(Math.cosh(g * x[0]), -2) - u, -1.0},
						{1.0, -(k + u)}
				})*/
		};
	}
	
	// method for evaluating Df(x, k) for bifurcation continunity
	private double [] [] Dfc (double [] x, double u)
	{
		double sech = Math.pow(Math.cosh(g * x[0]), -2.0);
		
		return new double [][] {
				{g * sech, -1.0, 0.0},
				{1, x[2], -x[1]},
				{2 * g * g * (x[2] + u) * Math.tanh(g * x[0]) * sech, 0, u - g * sech}
		};
	}
	
	// shooting method
	private double[] shootForFixedPoint (double [] v, double k, double error)
	{
		// local variable
		double [] v0 = new double [v.length];
		for (int i = 0; i < v.length; ++i)
			v0[i] = v[i];
		
		double [] fx = new double [v.length];
		double [][] Dfx;
		int counter = 0;
		
		try {
			// perform the iteration
			while (++counter < 100) {
				// evaluate f(x), Df(x) and new v
				fx = f(v0, k);
				Dfx = Df(v0, k);
				fx = Matrix.subtract(v0, Matrix.product(Matrix.inverse(Dfx), fx));
				
				// evaluate the percentage error in v1
				if (Matrix.getMaxError(fx, v0) < error)
					return fx;
				else
					v0 = fx;
			}
			throw new MatrixException ("Failure to converge");
		} catch (MatrixException e) {
			// TODO Auto-generated catch block
			System.out.println(e.getMessage());
		}
		return v0;
	}
	
	// method to evaluate the f(x)
	private double [] f(double [] x, double k)
	{
		return new double [] {
				-x[1] + Math.tanh(g * x[0]),
				x[0] - k * x[1]
		};
	}
	
	// method to evaluate the Jacobian of f(x)
	private double [] [] Df(double [] x, double k)
	{
		return new double [] [] {
				{g * Math.pow(Math.cosh(g * x[0]), -2), -1.0},
				{1.0, -k}
		};
	}
	
	/**
	 * This method displays a column matrix
	 * @param x
	 */
	private void display (double [] x)
	{
		for (int i = 0; i < x.length; ++i)
			System.out.println("\t|" + todp(x[i], 4) + "|");
	}
	
	/**
	 * This method returns a double value to a certain number of decimal places
	 * @param x The decimal value
	 * @param dp Decimal places
	 * @return
	 */
	private double todp (double x, int dp)
	{
		int value = (int) (x * Math.pow(10, dp) + 0.5);
		return (value / Math.pow(10, dp));
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new ContinunityTest ();
	}

}
