package phDModels;

import java.io.IOException;

import math.tools.Algebra;
import math.tools.Bifurcation;
import math.tools.ComplexNumber;
import math.tools.ExpNode;
import math.tools.Matrix;
import math.tools.MatrixException;

public class BVPLimitCycleStability {
	
	/**
	 * This method is basically used to analyze the limit
	 * cycle or the periodic stability of the system while computing the
	 * bifurcation point.
	 */
	public BVPLimitCycleStability () {
		// BVP Model
		/*
		StringBuffer [] model = new StringBuffer [] {
				new StringBuffer ("-y + tanh (g * x)"), 
				new StringBuffer("x - k * y")};
		String [] variables = new String [] {"g", "k", "x", "y"};
		String [] x = new String [] {"x", "y"};
		double [] values = new double [] {1.58, 0.841, 2.5, 0.0};
		//double [] values = new double [] {1.507, 0.855838, 2.5, 0.0};
		//double [] values = new double [] {1.521, 0.853, 2.5, 0.0};
		//double [] values = new double [] {1.034, 0.95, 2.5, 0.0};
		//double [] values = new double [] {2.0, 0.77, 2.5, 0.0};
		
		ExpNode [] f = Matrix.infixToExpNode(model, variables);
		
		// call the method to continuously analyze the limit cycles
		Bifurcation.showLimitCycleStability(null, f, x, variables, values, 
				"k", 0.5, 1.0, 1, 0.00001);
		*/
		
		// Ueta2004
		StringBuffer [] model = new StringBuffer [] {
				new StringBuffer ("-y1 + tanh (g * x1) - d * (x1 - x2)"), 
				new StringBuffer ("x1 - k1 * y1"),
				new StringBuffer ("-y2 + tanh (g * x2) - d * (x2 - x1)"), 
				new StringBuffer ("x2 - k2 * y2")};
		String [] variables = new String [] {
				"g", "k1", "k2", "d", "x1", "y1", "x2", "y2"};
		
		ExpNode [] f = Matrix.infixToExpNode(model, variables);
		
		Bifurcation.showLimitCycleStability(null, f, 
				new String [] {"x1", "y1", "x2", "y2"}, 
				new String [] {"g", "k1", "k2", "d", "x1", "y1", "x2", "y2"}, 
				new double [] {1.637, 1.0, 1.047, 0.337, 1.2, 1.5, 2.3, 1.8},
				"k2", 0.5, 0.8, 1, 0.001);
	}
	
	/**
	 * This method analyzes a range of k-values for the bifurcation point
	 * of the limit cycle stability and prints to the give file all the 
	 * computed stability for this range. The iteration process starts from the value given in 
	 * values and loops until a bifurcation point is found or an equilibrium is 
	 * found.
	 * 
	 * @param file The file into which the stability values are printed
	 * @param f The model equation
	 * @param x The independent variables of the state equation
	 * @param variables The variables in the model
	 * @param values The corresponding value of the variables
	 * @param vert The variable on the vertical axis
	 * @param hor The variable on the horizontal axis
	 * @param vert0 The starting point for the vertical axis
	 * @param vertf The end point for the vertical axis
	 * @param sectionIndex The Poincare section index in x
	 * @param dk k-value interval
	 * 
	 * @throws IllegalStateException If an equilibrium point is found.
	 */
	public void analyzeLimitCycle (java.io.File file, ExpNode [] f, String [] x, 
			String [] variables, double [] values, String vert, String hor, 
			double vert0, double vertf, int sectionIndex, double dk)
	{
		// initialize the number of the independent variables
		int n = x.length;
		double dt = 1E-2;
		
		// the list that has all the computed stability
		ComplexNumber [] u = new ComplexNumber [2 * n]; // also the variables
		
		// set the file and the writer
		java.io.FileWriter writer = setFileWriter (null);
					
		// set the file details
		try {
			// set the header
			writer.write("Limit Cycle Stability at g = " + 
					values [ExpNode.getIndex("g", variables)]);
			writer.write("\r\nk\t\th1\t\t\th2\t\t\tx\t\t\tt");
			writer.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}	
			
		// first define the time index and incorporate its variable name
		// and value
		int tIndex = variables.length;
		String t = "time";
		String [] var = new String [tIndex + 1];
		double [] val = new double [tIndex + 1];
		for (int i = 0; i < tIndex; ++i) {
			var [i] = variables[i];
			val [i] = values[i];
		}
		var [tIndex] = t;
		val [tIndex] = 0.0;
		
		// redefine variables to suit second variation equation
		// local variables
		var = Matrix.getVariationVar(var, x);
		double [] val1 = Matrix.getFirstVariationInitVal(val, x);
		ExpNode [] fx = Matrix.getVariationEq(f, x);
		
		// put the embedded variables into fx
		int [] [] Dhx = Matrix.Dhx(x, sectionIndex);
		int [] [] Dhu = Matrix.Dhu(x, sectionIndex);
		int [] Dqx = Matrix.Dqx (x, sectionIndex);
		
		// form the derivative function for the Newton-Raphson iteration
		ExpNode [] [] dw = Matrix.getFirstVariationXNode(x);
		ExpNode [] [] DFx = new ExpNode [n] [n];
		
		// first row element
		ExpNode [] [] DTu = Matrix.product(Dhx, Matrix.product(dw, Dhu));
		ExpNode [] DTt = Matrix.product(Dhx, f);;
		
		// complete DFx for the first row
		for (int i = 0; i < n - 1; ++i) {
			for (int j = 0; j < n - 1; ++j)
				// subtract the identity element from diagonal elements of DTu
				DFx [i][j] = (i == j ? ExpNode.makeBinOpNode ('-', DTu [i][j], 
						ExpNode.makeConstantNode(1.0)) : 
							DTu [i][j]);
			DFx [i][n - 1] = DTt[i];
		}
		
		// second row element
		ExpNode [] Dqu = Matrix.product(Dqx, Matrix.product(dw, Dhu));
		ExpNode Dqt = Matrix.product(Dqx, f);
		
		// complete second row
		for (int i = 0; i < n - 1; ++i)
			DFx [n - 1][i] = Dqu[i];
		DFx [n - 1] [n - 1] = Dqt;
		
		// Newton-raphson iteration
		double error = 1.0;
		String [] h = Matrix.hx(x, sectionIndex);
		String [] Dx = new String [n];
		int ttIndex = n - 1;
		for (int i = 0; i < h.length; ++i) 
			Dx[i] = h[i];
		Dx[ttIndex] = t;
		
		int [] DxIndices = ExpNode.getXIndex(Dx, var);
		int [] xIndices = ExpNode.getXIndex(
				Matrix.getVariationIndependentVar(x), var);
		double [] xx = new double [n];
		for (int i = 0; i < xx.length; ++i)
			xx [i] = val[DxIndices[i]];
		double [] ffx = new double [n];
		double [] [] Dff = new double [n][n];
		
		// start the analysis
		try {
			// loop through all the k-values
			for (double k = values [ExpNode.getIndex(vert, variables)];
					k >= vert0 && k <= vertf; k += dk) {
				val[ExpNode.getIndex(vert, var)] = k;
				
				// reset the correction iteration error
				error = 1.0;
				while (error > 1E-9) {
					try {
						// initialize the variables of the variational equation
						val1 = Matrix.getVariationInitVal (val, x);
						
						// get the next poincare section crossing
						val1 = Algebra.getNextCrossing(fx, xIndices, 
								sectionIndex, var, val1, tIndex, dt, 1.0E-12);
						// initialize xx
						for (int i = 0; i < n; ++i)
							xx[i] = val[DxIndices[i]];
						
						// reset ffx
						for (int i = 0; i < n - 1; ++i)
							ffx [i] = val1 [DxIndices[i]] - val[DxIndices[i]];
						ffx[n - 1] = val1[xIndices[sectionIndex]];
						
						Dff = Matrix.values(DFx, var, val1);
						
						ffx = Matrix.LUSolver(Dff, 
								Matrix.subtract(Matrix.product(Dff, xx), ffx));
						error = Matrix.getMaxError(ffx, xx);
						
						// reset the elements of Dx in val
						for (int i = 0; i < Dx.length; ++i)
							val[DxIndices[i]] = ffx[i];
						
					} catch (MatrixException e) {
						break;
					}
				}	// end while block
				
				// now evaluate the stability of the limit cycle
				ComplexNumber [] roots = Matrix.getEigenValues (
						Matrix.values(dw, var, val1));
				for (int l = 0; l < n; ++l)
					u [l] = roots[l];
				
				// put the Poincare Section crossing and the period
				for (int l = 0; l < DxIndices.length; ++l)
					u[n + l] = new ComplexNumber (val[DxIndices[l]]);
				
				// write the content of u to file
				writeToFile(writer, u, var, val1, dk);
			}
		} catch (IllegalStateException e) {
			
		}
	}
	
	/**
	 * This method tests the k-values of the given g-value for limit cycle
	 * bifurcation point. The iteration process starts from the value given in 
	 * values and loops until a bifurcation point is found or an equilibrium is 
	 * found.
	 * 
	 * @param fx
	 * @param var
	 * @param val
	 * @return
	 * @throws
	 */
	public double getBifurcationPoint (ExpNode fx, String [] var, double [] val)
	{
		
		return 0.0;
	}
	
	/**
	 * This method sets the file into which the limit cycle stability are
	 * printed and returns the writer
	 * @param file
	 */
	private static java.io.FileWriter setFileWriter (java.io.File file)
	{
		try {
			if (file == null) {
				int count = 0;
				file = new java.io.File ("C:\\Users\\Adesola\\Desktop\\limitBif.txt");
				while (file.exists()) {
					file = new java.io.File (
							"C:\\Users\\Adesola\\Desktop\\limitBif" + ++count + ".txt");
				}
				file.createNewFile();
			}
			
			return new java.io.FileWriter(file);
		} catch (IOException e) {
			e.printStackTrace();
			return null;
		}
	}
	
	/**
	 * This method prints the bifurcation parameters to file.
	 * 
	 * @param writer The java.io.FileWriter
	 * @param eigen The limit cycle stability values
	 * @param var The variables of the state model
	 * @param val The corresponding values of the state variables
	 * @param dk The incremental vertical value
	 */
	private static void writeToFile (java.io.FileWriter writer, 
			ComplexNumber [] eigen,	String [] var, double [] val, double dk) 
	{
		// begin writing into the file
		try {
			// print result to file
			writer.write(String.format ("\r\n%f", 
					val[ExpNode.getIndex("k", var)]));
			for (int j = 0; j < eigen.length; ++j)
				writer.write(String.format("\t%s", eigen[j]));
			writer.flush();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new BVPLimitCycleStability ();
	}

}
