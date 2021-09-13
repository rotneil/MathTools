package bifurcation;

import java.io.IOException;

import math.tools.Matrix;
import math.tools.MatrixException;

/**
 * This class is used to plot the complete two-parameter bifurcation diagram
 * for the BVP oscillator using the already discovered fixed points and their 
 * corresponding eigen values from BVPFixedPointBifurcation.
 * 
 * 
 * @author Nehemiah Oluwafemi
 *
 */
public class BVPBifurcation {
	// instance variable
	private double [] [] bifurcation = new double [2001][5];
	
	// no-argument constructor
	public BVPBifurcation() {
		/*
		 * This line of codes is used to confirm the eigenvalues at the pitchfork
		 * bifurcation points of bvp oscillator. It is useful for computing the
		 * the behaviour of a certain fixed point over a range of g-values
		 *//*
		try {
			double g = 1.8, sech, error;
			double [] x = new double [2];
			double [] xn = new double [2];
			double [] f = new double [2];
			double [] [] Df = new double [2][2];
			
			// choose file location
			java.io.File file = new java.io.File ("C:\\Users\\Adesola\\Desktop\\" +
					"bifAnly1.txt");
			int count = 1;
			while (file.exists()) {
				file = new java.io.File (
						"C:\\Users\\Adesola\\Desktop\\bifAnly" + ++count + ".txt");
			}
			file.createNewFile();
			java.io.FileWriter writer = new java.io.FileWriter(file);
			
			// prepare the heading
			writer.write("Analysis at g = " + g + "\r\n");
			
			// evaluate the stability of the fixed points
			for (double k = 1.0 / g; k <= 1.0; k += 0.0001) {
				error = 1.0;
				//x = new double [] {-1.5, 0.8};		// for origin
				x = new double [] {1.2133, 1.3867};	// for a certain fixed point
				
				while (error > 1.0E-6) {
					sech = 1.0 / Math.cosh(g * x[0]);
					f = new double [] {
							-x[1] + Math.tanh(g * x[0]),
							x[0] - k * x[1]};
					Df = new double [] [] {
							{g * sq(sech), -1.0},
							{1.0, -k}};
					xn = Matrix.subtract(x, Matrix.product(Matrix.inverse(Df), f));
					error = Matrix.getMaxError(xn, x);
					x = xn;
				}
				ComplexNumber [] h = Matrix.getEigenValues(Df);
				
				writer.write(String.format("x = %f4\ty = %f4\tk = %f6\t" +
						"h1 = %f6 " + (h[0].getImaginaryNumber() >= 0.0? "+ i":"- i") +
						"%f6\tmod = %f6\r\n", x[0], x[1], k, h[0].getRealNumber(), 
						h[0].getImaginaryNumber(), h[0].getModulus()));
				writer.flush();
			}

			writer.close();
		} catch (java.io.IOException e) {
			e.printStackTrace();
		}
		catch (MatrixException e) {
			e.printStackTrace();
		}*/
		
		// fill in the continuity patterns
		double dg = 0.001;
		getAndronovCont (0, 0.0, 0.0, 0.8, 0.8, 0.6, dg);/*
		getAndronovCont (0, 0.0, 0.0, 0.82, 0.82, 0.572364, -dg);
		getPitchForkCont (1, 0.0, 0.0, 1.22, 0.82, 0.0, dg);
		/*getPitchForkCont (1, 0.0, 0.0, 1.22, 0.82, 0.0, -dg);
		double [] x = refinePoint(0.56549, 0.769003, 1.8, 0.735356, 0.00009, 0.67758);
		getBifurcation(2, x[0], x[1], 1.8, x[2], x[4], dg);
		getBifurcation(2, x[0], x[1], 1.8, x[2], x[4], -dg);
		// use the bifurcation points from BVPLimitCycleBifurcation.java to 
		// compute the tangent continuity bifurcation for the limit cycles
		getTangentCont (3, 0.454580693582, 1.8, 0.80375950, 17.094267423629763, dg);
		getTangentCont (3, 0.454580693582, 1.8, 0.80375950, 17.094267423629763, dg);
		getTangentCont (3, 0.287080, 1.522, 0.852958, 20.197709, -dg);/*
		getTangentCont (3, 0.740276, 1.4, 0.878286, 20.026798, dg);
		getTangentCont (3, 0.740276, 1.4, 0.878286, 20.026798, -dg);/*
		getTangentCont (4, 0.734205, 1.8, 0.786556, 12.966309, dg);
		getTangentCont (4, 0.734205, 1.8, 0.786556, 12.966309, -dg);*/
		
		// print bifurcation points to file
		printBifurcation ();
	}
	
	// this method computes the continuity params for tangent bifurcation of
	// the limit cycle
	private void getTangentCont (int index, double x0, double g, double k,
			double t, double dg)
	{
		// intialize value
		double [] x = new double [] {x0, t, k};
		double [] xx = new double [19];
		double [] f = new double [x.length];
		double [] [] Df = new double [x.length] [x.length];
		double f1, f2;
		double dt = 0.001;
		double error = 1.0;
		
		// do the iteration
		try {
			for (double gg = g; gg > 1.0 && gg <= 2.0; gg += dg) {
				try {
					error = 1.0;
					while (error > 1.0E-8) {
						// dynamics over a cycle at the poincare section crossing
						xx = getCycleDynamics (x[0], gg, x[2], dt);
						x[1] = xx[18];		// adjust the time to poincare crossing
						
						// evaluate params for Newton's interpolation
						f1 = -xx[1] + Math.tanh(gg * xx[0]);
						f2 = xx[0] - x[2] * xx[1];
						f =  new double [] {xx[0] - x[0], xx[1], 2 - xx[2] - xx[5]};
						Df = new double [] [] {
							{xx[2] - 1.0, f1, xx[6]},
							{xx[3], f2, xx[7]},
							{-xx[8] - xx[11], 
								-f1 * (xx[8] + xx[11]) - f2 * (xx[10] + xx[13]), 
								-xx[6] * (xx[8] + xx[11]) - xx[7] * (xx[10] + xx[13])}};
						
						// evaluate the current x, error and invert
						f = Matrix.subtract(x, Matrix.product(Matrix.inverse(Df), f));
						error = Matrix.getMaxError(f, x); //Math.abs(f[2] - x[2]);
						x = f;
					}
					System.out.printf("g = %f6\tx = %f6\tt = %f6\tk = %f6\n", gg, 
							x[0], x[1], x[2]);
					
					// assign the k-value at the specified index
					bifurcation[((int) (gg * 1000 + 0.5))] [index] = x[2];
				} catch (IllegalStateException e) {
					System.out.printf("g = %f6\t Jumping iterate\n", gg);
					continue;
				}
			}
		} catch (MatrixException e) {
			System.out.println ("Tangent Forward Continuity for limit cycle: k" + 
					(index + 1) + " completed.");
		}
	}
	
	/**
	 * This method is used to compute an iteration over a cycle of the system orbit.
	 * @param x0
	 * @param g
	 * @param k
	 * @return xx, which contains the dynamics of the system, i.e. first and second
	 * variational matrices of the system
	 */
	private double [] getCycleDynamics (double x0, double g, double k, double dt)
	{
		// local variables
		double [] v = new double [] {x0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
				0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		double [] vv = rungeKutta (v, g, k, dt);
		int count = 1;
		while (!(v[1] < 0.0 && vv[1] >= 0.0)){
			// check that orbit has not entered equilibrium point
			if (Math.abs(v[0] - vv[0]) <= 1.0E-12) {
				count = 1;
				v = new double [] {x0 + 0.01, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
						0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
				vv = rungeKutta (v, g, k, dt);
				//throw new IllegalStateException ("Orbit at equilibrium point");
			}
			v = vv;
			vv = rungeKutta (v, g, k, dt);
			++count;
		}
		
		return timeHalf (v, g, k, (count - 1) * dt, dt, 1.0E-18);
	}
	
	// timeHalf method for Poincare section mapping
	private double [] timeHalf (double [] v, double g, double k, double t, double dt,
			double err)
	{
		// local variable
		double tt = 0.0;
		double dtt = dt * 0.5;
		double [] vv = rungeKutta (v, g, k, dtt);
		while (Math.abs(vv[1]) > err) {
			if (vv[1] > 0.0)
				dt = dtt;
			else if (vv[1] < 0.0)
				tt = dtt;
			dtt = (dt + tt) * 0.5;
			vv = rungeKutta (v, g, k, dtt);
		}
		
		vv[18] = t + dtt;
		return vv;
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
 		double P = g * Math.pow((Math.cosh(g * v[0])), -2.0);
 		double Q = g * Math.tanh(g * v[0]);
 		
 		vv[0] = -v[1] + Math.tanh(g * v[0]);
 		vv[1] = v[0] - k * v[1];
 		
 		vv[2] = P * v[2] - v[3];
 		vv[3] = v[2] - k * v[3];
 		vv[4] = P * v[4] - v[5];
 		vv[5] = v[4] - k * v[5];
 		
 		vv[6] = P * v[6] - v[7];
 		vv[7] = v[6] - k * v[7] - v[1];
 		
 		vv[8] = P * v[8] - 2 * P * Q * sq(v[2]) - v[9];
 		vv[9] = v[8] - k * v[9];
 		vv[10] = P * v[10] - v[11] - 2 * P * Q * v[2] * v[4];
 		vv[11] = v[10] - k * v[11];
 		vv[12] = P * v[12] - v[13] - 2 * P * Q * sq(v[4]);
 		vv[13] = v[12] - k * v[13];
 		
 		vv[14] = P * v[14] - v[15] - 2 * P * Q * v[2] * v[6];
 		vv[15] = v[14] - k * v[15] - v[3];
 		vv[16] = P * v[16] - v[17] - 2 * P * Q * v[4] * v[6];
 		vv[17] = v[16] - k * v[17] - v[5];
 		
 		return vv;
 	}
	
 	// method to return the continuity for Andronov hopf bifurcation
	private void getAndronovCont (int index, double x1, double x2, double g, 
			double k, double im, double dg)
	{
		// intialize value
		double [] x = new double [] {x1, x2, k, im};
		double [] f = new double [x.length];
		double [] [] Df = new double [x.length] [x.length];
		double sech = 0.0;
		double error = 1.0;
		
		// do the iteration
		try {
			for (double gg = g; gg >= 0.0 && gg <= 2.0; gg += dg) {
				error = 1.0;
				while (error > 1.0E-6) {
					// evaluate differentials
					sech = 1.0 / Math.cosh (gg * x[0]);
					f =  new double [] {
							-x[1] + Math.tanh(g * x[0]),
							x[0] - x[2] * x[1],
							-x[3] * x[3] - gg * x[2] * sech * sech + 1,
							x[3] * (x[2] - gg * sech * sech)};
					Df = new double [] [] {
							{gg * sech * sech, -1.0, 0.0, 0.0},
							{1.0, -x[2], -x[1], 0.0},
							{2 * gg * gg * x[2] * Math.tanh(gg * x[0]) * sech * sech, 0, 
								gg * sech * sech, -2.0 * x[3]},
							{2 * x[3] * gg * gg * Math.tanh(gg * x[0]) * sech * sech, 0, x[3],
									x[2] - gg * sech * sech
							}};
					
					// evaluate the current x, error and invert
					f = Matrix.subtract(x, Matrix.product(Matrix.inverse(Df), f));
					error = Matrix.getMaxError(f, x);
					x = f;
				}
				// assign the k-value at the specified index
				bifurcation[((int) (gg * 1000 + 0.5))] [index] = x[2];
			}
		} catch (MatrixException e) {
			System.out.println ("Andronov Forward Continuity for k" + 
					(index + 1) + " completed.");
		}
	}
	
	// method to refine the bifurcation point
	private double [] refinePoint (double x1, double x2, double g, 
			double k, double r, double m)
	{
		double [] x = new double [] {x1, x2, k, r, m};
		double [] f = new double [x.length];
		double [] [] Df = new double [x.length] [x.length];
		double sech = 0.0, tanh = 0.0;
		
		// do the iteration
		try {
			while (Math.abs(x[3]) > 1.0E-7) {
				// evaluate differentials
				sech = 1.0 / Math.cosh (g * x[0]);
				tanh = Math.tanh(g * x[0]);
				f =  new double [] {
						-x[1] + tanh,
						x[0] - x[2] * x[1],
						1 - (x[2] + x[3]) * (g * sq(sech) - x[3]) - sq(x[4]),
						x[4] * (x[2] + 2 * x[3] - g * sq(sech)),
						x[3]};
				
				Df = new double [] [] {
						{g * sq(sech), -1.0, 0.0, 0.0, 0.0},
						{1.0, -x[2], -x[1], 0.0, 0.0},
						{2 * sq(g) * (x[2] + x[3]) * tanh * sq(sech), 0, 
						g * sq(sech) - x[3], x[2] + 2 * x[3] - g * sq(sech), -2 * x[4]},
						{2 * x[4] * sq(g) * tanh * sq(sech), 0.0, x[4], 2 * x[4],
							x[2] + 2 * x[3] - g * sq(sech)},
						{0.0, 0.0, 0.0, 1.0, 0.0}};
				
				// evaluate the current x, error and invert
				x = Matrix.subtract(x, Matrix.product(Matrix.inverse(Df), f));
			}
		} catch (MatrixException e) {
			System.out.println();
		}
		return x;
	}
	
	// method to return the continuity for Andronov hopf bifurcation
	private void getBifurcation (int index, double x1, double x2, double g, 
			double k, double m, double dg)
	{
		// intialize value
		double [] x = new double [] {x1, x2, k, m};
		double [] f = new double [x.length];
		double [] [] Df = new double [x.length] [x.length];
		double sech = 0.0, tanh = 0.0;
		double error = 1.0;
		
		// do the iteration
		try {
			// start the continuity equation
			for (double gg = g; gg >= 1.0 && gg <= 2.0; gg += dg) {
				error = 1.0;
				while (error > 1.0E-6) {
					// evaluate differentials
					sech = 1.0 / Math.cosh (gg * x[0]);
					tanh = Math.tanh(gg * x[0]);
					f =  new double [] {
							-x[1] + tanh,
							x[0] - x[2] * x[1],
							1.0 - x[2] * gg * sq(sech) - sq(x[3]),
							x[3] * (x[2] - gg * sq(sech))};
					Df = new double [] [] {
							{gg * sq(sech), -1.0, 0, 0},
							{1.0, -x[2], -x[1], 0.0},
							{2 * sq(gg) * x[2] * tanh * sq(sech), 0, -gg * sq(sech), -2 * x[3]},
							{2 * x[3] * sq(gg) * tanh * sq(sech), 0, x[3], x[2] - gg * sq(sech)}};
					
					// evaluate the current x, error and invert
					f = Matrix.subtract(x, Matrix.product(Matrix.inverse(Df), f));
					error = Matrix.getMaxError(f, x);
					x = f;
				}
				// assign the k-value at the specified index
				bifurcation[((int) (gg * 1000 + 0.5))] [index] = x[2];
			}
		} catch (MatrixException e) {
			System.out.println ("Andronov Forward Continuity for k" + 
					(index + 1) + " completed.");
		}
	}
	
	// method to return the conitinuity pattern for Pitchfork bifurcation
	private void getPitchForkCont (int index, double x1, double x2, double g, 
			double k, double u, double dg)
	{
		// intialize value
		double [] x = new double [] {x1, x2, k};
		double [] f = new double [x.length];
		double [] [] Df = new double [x.length] [x.length];
		double sech = 0.0;
		double error = 1.0;
		
		// do the forward iteration
		for (double gg = g; gg >= 0.0 && gg <= 2.0; gg += dg) {
			error = 1.0;
			while (error > 1.0E-6 && x[2] >= 0 && x[2] <= 2.0) {
				// evaluate differentials
				sech = 1.0 / Math.cosh (gg * x[0]);
				f =  new double [] {
						-x[1] + Math.tanh(gg * x[0]),
						x[0] - x[2] * x[1], 
						1.0 - (x[2] + u) * (gg * sech * sech - u)};
				Df = new double [][] {
						{gg * sech * sech, -1.0, 0.0},
						{1, -x[2], -x[1]},
						{2 * gg * gg * (x[2] + u) * Math.tanh(gg * x[0]) 
							* sech * sech, 0.0, u - gg * sech * sech}};
				
				// evaluate the current x, error and invert
				try {
					f = Matrix.subtract(x, Matrix.product(Matrix.inverse(Df), f));
					error = Matrix.getMaxError(f, x);
					x = f;
				} catch (MatrixException e) {
					break;
				}
			}
			// assign the k-value at the specified index
			bifurcation[((int) (gg * 1000 + 0.5))] [index] = x[2];
		}
		System.out.println ("Pitchfork " + (dg > 0.0 ? "Forward" : "Backward") +
				" Continuity for k" + (index + 1) + " completed.");
	}
	
	// method to write points to file
	public void printBifurcation ()
	{
		try {
			// choose file location
			java.io.File file = new java.io.File ("C:\\Users\\Adesola\\Desktop\\bif1.txt");
			int count = 1;
			while (file.exists()) {
				file = new java.io.File (
						"C:\\Users\\Adesola\\Desktop\\bif" + ++count + ".txt");
			}
			file.createNewFile();
			java.io.FileWriter writer = new java.io.FileWriter(file);
			
			// prepare the heading
			writer.write("g\tK1\tK2\tK3\tK4\tK5\r\n");
			
			// output values to file
			for (int i = 0; i < bifurcation.length; ++i)
				writer.write(String.format("%f4\t%f6\t%f6\t%f6\t%f6\t%f6\r\n",
						(i * 0.001), bifurcation[i][0], bifurcation[i][1], 
						bifurcation[i][2], bifurcation[i][3], bifurcation[i][4]));
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private double sq (double x) { return x * x; }
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new BVPBifurcation();
	}

}
