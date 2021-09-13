package bifurcation;

import java.io.IOException;

import math.tools.ComplexNumber;
import math.tools.Matrix;
import math.tools.MatrixException;

public class BVPLimitCycleBifurcation {

	public BVPLimitCycleBifurcation() {
		// TODO Auto-generated constructor stub
		try {
			// local variables
			double g = 1.3;
			double dt = 0.001;
			
			// choose file location
			java.io.File file = new java.io.File ("C:\\Users\\Adesola\\Desktop\\" +
					"lcBif1.txt");
			int count = 1;
			while (file.exists()) {
				file = new java.io.File (
						"C:\\Users\\Adesola\\Desktop\\lcBif" + ++count + ".txt");
			}
			file.createNewFile();
			java.io.FileWriter writer = new java.io.FileWriter(file);
			
			// prepare the heading
			writer.write("Analysis at g = " + g + "\r\n");
			
			// iterate through values of k
			double [] x = new double [7];
			double [] xn = new double [2];
			double [] f = new double [2];
			double [] [] Df = new double [2] [2];
			double error = 1.0;
			for (double k = 0.87890544; k < 1.0; k += 1.0E-15) {
				// iterate variables
				xn = new double [] {4.0};
				
				// shoot for fixed point
				error = 1.0;
				while (error > 1.0E-12) {
					// evaluate variables
					x = getNextCrossing (xn[0], g, k, dt);
					f = new double [] {x[0] - xn[0], x[1]};
					Df = new double [] [] {
							{x[2] - 1.0, -x[1] + Math.tanh(g * x[0])},
							{x[3], x[0] - k * x[1]}};
					
					// evaluate the new point
					f = Matrix.subtract(xn, Matrix.product(Matrix.inverse(Df), f));
					error = Matrix.getMaxError(f, xn);
					xn = f;
				}	// end while
				
				// evaluate the chi at the fixed point
				ComplexNumber [] u = Matrix.getEigenValues(
						new double [][] {{x[2], x[4]}, {x[3], x[5]}});
				
				writer.write(String.format("x = %s\tk = %s\tt = %s\t" +
						"h1 = %f4 " + (u[0].getImaginaryNumber() >= 0.0? "+ i":"- i") +
						"%f4\t" +
						"h2 = %f4 " + (u[1].getImaginaryNumber() >= 0.0? "+ i":"- i") +
						"%f4\r\n", x[0], k, x[6], u[0].getRealNumber(), 
						u[0].getImaginaryNumber(), u[1].getRealNumber(), 
						u[1].getImaginaryNumber()));
				writer.flush();
			}	// end for (k)
			writer.close();
		} catch (MatrixException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (IllegalStateException e) {
			
		}
	}
	
	/* method to get the next crossing of the Poincare section for the x-value
	** initially on the section. The return value has the time required to 
	** make this crossing
	*/
	private double [] getNextCrossing (double x, double g, double k, double dt)
	{
		// local variables
		double [] v = new double [] {x, 0.0, 1.0, 0.0, 0.0, 1.0};
		double [] vv = rungeKutta (v, g, k, dt);
		int count = 1;
		while (!(v[1] < 0.0 && vv[1] >= 0.0)){
			// check that the orbit is not at equilibrium point
			if (Math.abs(v[0] - vv[0]) < 1.0E-12)
				throw new IllegalStateException ("The orbit is at equilibrium point.");
			v = vv;
			vv = rungeKutta (v, g, k, dt);
			
			++count;
		}
		
		/*// interpolate the values of points to move them to Poincare section
		double a1 = v[1], a2 = vv[1];
		double [] p = new double [v.length + 1];
		
		for (int i = 0; i < v.length; ++i)
			p[i] = a2 * v[i] / (a2 - a1) + a1 * vv[i] / (a1 - a2);
		
		// evaluate the time into p
		p[v.length] = (count - a2 / (a2 - a1)) * dt;
		
		return p;
		*/
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
		
		// initialize the return variable
		double [] p = new double [v.length + 1];
		
		for (int i = 0; i < v.length; ++i)
			p[i] = vv[i];
		
		// evaluate the time into p
		p[v.length] = t + dtt;
		
		return p;
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
 		double p = g * Math.pow(Math.cosh(g * v[0]), -2.0);
 		double q = g * Math.tanh(g * v[0]);
 		
 		vv[0] = -v[1] + Math.tanh(g * v[0]);
 		vv[1] = v[0] - k * v[1];
 		
 		vv[2] = p * v[2] - v[3];
 		vv[3] = v[2] - k * v[3];
 		vv[4] = p * v[4] - v[5];
 		vv[5] = v[4] - k * v[5];
 		
 		return vv;
 	}
 	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new BVPLimitCycleBifurcation();
	}

}
