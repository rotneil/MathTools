package phDModels;

import java.io.IOException;

import math.tools.Algebra;
import math.tools.ComplexNumber;
import math.tools.ExpNode;
import math.tools.Matrix;
import math.tools.MatrixException;

/**
 * This program is used to identify the bifurcation points of the equilibrium 
 * by analyzing the stability of the equilibrium points of the models of my project. 
 * 
 * The results will then be used in the two parameter bifurcations for the model
 * 
 * 
 * @author Oluwafemi Nehemiah
 *
 */
public class StabilityAnalysis extends javax.swing.JPanel
{
	private static final long serialVersionUID = 1L;
	double k1 = 0.82, k = 0.82,
			g2 = 0.62, k2 = 0.98, s = 1.0E-3, e1 = -5.0E-5, e3 = 0.005;
	
	private int modelNumber;
	StringBuffer [] modelEquation;
	String [] modelVariables, x;
	ExpNode [] fx;
	ExpNode [] [] Df;
	ComplexNumber [] h;
	ComplexNumber [] [] dir;
	
	java.io.FileWriter writer;
	java.io.File file;
	
	// constructor
	public StabilityAnalysis() {
		// TODO Auto-generated constructor stub
		// iterate through the vertical i.e values of d
		modelNumber = 1;
		modelEquation = getModelEquation (modelNumber);
		modelVariables = getModelVariables (modelNumber);
		x = getModelIndependentVar (modelNumber);
		
		fx = Matrix.infixToExpNode(modelEquation, modelVariables);
		Df = Matrix.jacobian(fx, x);
		
		h = new ComplexNumber [fx.length];
		dir = new ComplexNumber [2][fx.length];
		
		// instantiate file and its writer
		int count = 1; 
		try {
			file = new java.io.File ("C:\\Users\\Adesola\\Desktop\\" +
					"m" + modelNumber + "bif1.txt");
			while (file.exists()) {
				file = new java.io.File (
						"C:\\Users\\Adesola\\Desktop\\m" + modelNumber +
						"bif" + ++count + ".txt");
			}
			file.createNewFile();
			writer = new java.io.FileWriter(file);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		model1Analysis ();
		//model2Analysis ();
	}
	
	private void model2Analysis ()
	{
		// the analysis of the system is done along the diagonal of axes k1-k2
		// plane
		try {
			writer.write("Model " + modelNumber + " Analysis " +
					"of the equilibrium at \r\ng = 1.637 and d = 0.337 " + "\r\n" +
					"Intended for k1 against k2 bifurcation plot\r\n");
			writer.flush();			
			
			// begin the equilibrium bifurcation test
			// evaluate the jacobian of the system
			ExpNode [] [] Df = Matrix.jacobian(fx, x);
			
			// iterate through different values of delta
			for (double k1 = 0.0, k2 = 0.0; k1 <= 2.0; k1 += 0.0001) {
				k2 = k1;
				h = Matrix.getEigenValues(Matrix.values(Df, modelVariables, 
					new double [] {0.0, 0.0, 0.0, 0.0, 1.637, k1, k2, 0.337}));
				
				// check the eigen values but print the line with bifurcation
				ComplexNumber [] bif = compareEigenValues (dir, h);
				
				if (bif != null) {
					writer.write (String.format("\r\n%.4f\t%.4f\t\t",
							k1, k2));
					for (int ij = 0; ij < h.length; ++ij)
						writer.write(String.format("h%d = %s  ", 
								(1 + ij), str(h[ij])));
					writer.flush();
				}
			}
			
			// perform the limit cycle bifurcation test
			analyzeLimitCycle (fx, x, 3, "k2", "k1", 0.3, 1.8, 0.3, 1.8, 0.3, 1E-2);
			
		} catch (MatrixException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	// analysis of model 1
	private void model1Analysis ()
	{
		int count = 1; // number of different initial condition trials
		
		try {
			// prepare the heading
			writer.write("Model " + modelNumber + " Analysis at\r\ng2 = " + 
					g2 + "\r\nk1 = " + k1 +	"\r\nk2 = " + k2 + "\r\ns = " + s + 
					"\r\ne1 = " + e1 + "\r\ne3 = " + e3 + "\r\n" +
					"Intended for d against g1 bifurcation plot\r\n");
			writer.flush();
			
			// iterate through different values of delta
			//for (double d = 0.0; d <= 2.0; d += 0.05) {
			for (double k = 0; k <= 2.0; k += 0.05) {
				// write the heading
				writer.write("\r\n\r\nAnalysis at delta = " + k ); // d);
				writer.flush();
				
				// initialize x to the origin
				count = 0;
				// iterate through different values of g
				for (double g = 0.0; g <= 2.0; g += 0.001) {
					h = Matrix.getEigenValues(Matrix.values(Df, modelVariables, 
							//new double [] {0.0, 0.0, 0.0, 0.0, g, k, d}));
							new double [] {0.0, 0.0, g, k}));
					
					// check the eigen values but print the line with bifurcation
					ComplexNumber [] bif = compareEigenValues (dir, h);
					
					if (bif != null) {
						writer.write(String.format("\r\ng = %.3f  ", (g - 0.001)));
						for (int ij = 0; ij < bif.length; ++ij)
							writer.write(String.format("h%d = %s  ", 
									(1 + ij), str(bif[ij])));
						++count;
					}
				}
				writer.write("\r\nAnalysis Points: " + count + " Counts");
				
				// iterate through different values of g
				for (double g = 1.53; g <= 1.54; g += 0.000001) {
					for (int i = 0; i < modelEquation.length; ++i)
						modelEquation[i].append(" - x" + (i + 1));
					ExpNode [] fxx = Matrix.infixToExpNode(modelEquation, modelVariables);
					ExpNode [][] dfxx = Matrix.jacobian(fxx, x);
					double [] val = new double[] {1.5, 1.5, g, k};
					double [] xx = new double [] {1.5, 1.5};
					double [] ffxx = new double [fxx.length];
					double [] [] dffxx =new double [fxx.length][fxx.length];
					double error = 1.0;
					while (error >= 1.0E-6) {
						ffxx = Matrix.values(fxx, modelVariables, val);
						dffxx = Matrix.values(dfxx, modelVariables, val);
						ffxx = Matrix.subtract(xx,
								Matrix.product(Matrix.inverse(dffxx), ffxx));
						error = Matrix.getMaxError(ffxx, xx);
						xx = ffxx;
						val[0] = xx[0];
						val[1] = xx[1];
					}
					h = Matrix.getEigenValues(Matrix.values(Df, modelVariables, 
							//new double [] {0.0, 0.0, 0.0, 0.0, g, k, d}));
							new double [] {xx[0], xx[1], g, k}));
					
					// check the eigen values but print the line with bifurcation
					ComplexNumber [] bif = compareEigenValues (dir, h);
					
					if (bif != null) {
						writer.write(String.format("\r\ng = %.3f  x1 = %.4f  " +
								"x2 = %.4f  ", (g - 0.001), xx[0], xx[1]));
						for (int ij = 0; ij < bif.length; ++ij)
							writer.write(String.format("h%d = %s  ", 
									(1 + ij), str(bif[ij])));
						++count;
					}
				}
				writer.flush();
			}
			
			// perform limitcycle analysis tes
			analyzeLimitCycle (fx, x, 1, "k", "g", 0.6, 0.82, 1.699, 2.0,
					0.5, 0.01);
			
			writer.close();
		} catch (MatrixException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void analyzeLimitCycle (ExpNode [] fx, String [] x, int section,
			String yAxis, String xAxis, double y1, double y2, double x1, double x2,
			double dy, double dx) throws IOException, MatrixException
	{
		// get the independent variable indices in model variables
		ExpNode [] ffx = Matrix.getVariationEq(fx, x);
		String [] var = Matrix.getVariationVar (modelVariables, x);
		x = Matrix.getVariationIndependentVar(x);
		int [] indices = ExpNode.getXIndex(x, var);
		double dt = 0.01, err = 1.0E-18;
		
		writer.write("\r\n\r\nLimit Cycle Analysis\r\n");
		writer.write(String.format("At " + yAxis + " = %.6f\r\n", y1));
			writer.flush();
		
		for (double kk = y1; kk <= y2; kk += dy) {
			for (double gg = x1; gg <= x2; gg += dx) {
				double [] val = new double [] {
						4.0, 1.2, -0.5, 0.0, gg, kk, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
						0.0, 0.0, 1.0, 0.0,	0.0, 0.0, 0.0, 1.0, 0.0};
				double [] xx = new double [] {val[0], 1.0};
				double [] ff = new double [fx.length];
				double [] [] Dff = new double [fx.length][fx.length];
				double error = 1.0;
				
				while (error > 1.0E-10) {
					val = new double [] {val[0], val[1], val[2], 0.0, gg, kk,	
							1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
							0.0, 0.0, 1.0, 0.0,	0.0, 0.0, 0.0, 1.0, 0.0};
					val = Algebra.getNextCrossing(ffx, 
							indices, section, var, val, -1, dt, err);
					ff = new double [] {
							val[0] - xx[0], val[section] };
					Dff = new double [] [] {
							{val[4] - 1, ExpNode.getValue(fx[0], var, val)},
							{val[6], ExpNode.getValue(fx[1], var, val)}};
					
					// evaluate new point
					ff = Matrix.subtract(xx, 
							Matrix.product(Matrix.inverse(Dff), ff));
					error = Matrix.getMaxError(ff, xx);
					val[0] = xx[0] = ff[0];
					xx[1] = val[1];
				}
				h = Matrix.getEigenValues(new double [] [] {
						{val[4], val[5]}, {val[6], val[7]}});
				writer.write(String.format("\r\ng = %.6f   " +
						"x = %.6f  t = %.6f  ", gg, xx[0], 
						val[val.length - 1]));
				for (int ij = 0; ij < h.length; ++ij)
					writer.write(String.format("h%d = %s  ", 
							(1 + ij), str(h[ij])));
				writer.flush();		
			}
		}
	}
	
	private String str (ComplexNumber h)
	{
		return String.format(
				(h.getRealNumber() >= 0 ? " " : "-") + "%.6f" + 
				(h.getImaginaryNumber() >= 0 ? " + i" : " - i") + "%.6f",
				Math.abs(h.getRealNumber()), Math.abs(h.getImaginaryNumber()));
	}
	
	// method that compares the direction of the eigen values
	private ComplexNumber[] compareEigenValues (
			ComplexNumber[][] dir, ComplexNumber[] h)
	{
		// local variable
		boolean realBif = false;
		boolean imBif = false;
		
		// check against null elements
		if (dir [0][0] != null) {
			// iterate through the eigen values
			for (int i = 0; i < h.length; ++i) {
				// compare real
				if (dir [0][i].getRealNumber() * h[i].getRealNumber() < 0.0)
					realBif = true;
				// compare imaginary
				if (dir [0][i].getImaginaryNumber() * h[i].getImaginaryNumber() < 0.0)
					imBif = true;
			}
		}
		
		// update the complex direction
		dir [0] = dir [1];
		dir [1] = h;
		
		if (realBif || imBif)
			return dir [0];
		
		return null;
	}
	
	// method to return the model equation
	private StringBuffer [] getModelEquation (int modelName)
	{
		switch (modelName) {
		case 1:		// Model 17
			return new StringBuffer [] {
					new StringBuffer ("-x2 + tanh (g * x1)"),
					new StringBuffer ("x1 - k * x2")
			};
			
		case 2:		// Ueta2004
			return new StringBuffer [] {
					new StringBuffer ("-x2 + tanh (g * x1) - d * (x1 - x3)"),
					new StringBuffer ("x1 - k1 * x2"),
					new StringBuffer ("-x4 + tanh (g * x3) - d * (x3 - x1)"),
					new StringBuffer ("x3 - k2 * x4")
			};
		case 3:
			return new StringBuffer [] {
					new StringBuffer ("-x2 + tanh (g * x1)"),
					new StringBuffer ("x1 - k * x2 + d * k * (x2 - x4)"),
					new StringBuffer ("-x4 + tanh (g * x3)"),
					new StringBuffer ("x3 - k * x4 + d * k * (x4 - x2)")
			};
		case 4:
			return new StringBuffer [] {
				new StringBuffer ("-x2 + tanh (g1 * x1) - " +
						"d * (e1 + e3 * x5 ^ 2) * (x1 - x3)"),
				new StringBuffer ("x1 - k1 * x2"),
				new StringBuffer ("-x4 + tanh (g2 * x3) - " +
						"d * (e1 + e3 * x5 ^ 2) * (x3 - x1)"),
				new StringBuffer ("x3 - k2 * x4"),
				new StringBuffer ("s * (x1 - x3)")
			};
		default:
			return null;
		}
	}
	
	// method to return the variable for the model
	private String [] getModelVariables (int modelName)
	{
		switch (modelName) {
		case 1:
			return new String [] {"x1", "x2", "g", "k"};
		case 2:
			return new String [] {"x1", "x2", "x3", "x4", "g", "k1", "k2", "d"};
		case 3:
			return new String [] {"x1", "x2", "x3", "x4", "g", "k", "d"};
		case 4:
			return new String [] {"x1", "x2", "x3", "x4", "x5", "g1", "g2",
					"k1", "k2", "d", "s", "e1", "e3"};
		default:
			return null;
		}
	}
	
	// method to return the independent variables
	private String [] getModelIndependentVar (int modelName)
	{
		switch (modelName) {
		case 1: // model 17
			return new String [] {"x1", "x2", "x3", "x4", "x5"};
		case 2: case 3:
			return new String [] {"x1", "x2", "x3", "x4"};
		case 4:
			return new String [] {"x1", "x2"};
		default:
			return null;
		}
	}
	
	/**
	 * This method returns a random number
	 * @return
	 */
	private double random () { return -3 + Math.random() * 6; }
	private double randomPhi () { return -1 + Math.random() * 2; }
	private double sq (double x) { return x * x; }
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new StabilityAnalysis();
	}

}
