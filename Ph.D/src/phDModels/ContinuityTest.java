package phDModels;

import java.io.IOException;

import math.tools.Algebra;
import math.tools.Bifurcation;
import math.tools.ComplexNumber;
import math.tools.ExpNode;
import math.tools.Matrix;
import math.tools.MatrixException;

public class ContinuityTest {
	
	// instance variables
	private static final long serialVersionUID = 1L;
	java.io.FileWriter mWriter;
	//double k1 = 0.82, g2 = 0.62, k2 = 0.98, s = 1.0E-3, e1 = -5.0E-5, e3 = 0.005;
	double xmax = 2.0, xmin = 0.0;
	double ymax = 2.0, ymin = 0.0;
	
	private int modelNumber = 1;
	private double [] [] kValue;
	private double dg = 0.001;
	
	private StringBuffer [] equation;
	private String [] var;
	private String [] x;
	private ExpNode [] f, fx;
	private ExpNode [] [] Dfx;
	
	public ContinuityTest() {
		// TODO Auto-generated constructor stub
		kValue = new double [(int)(2.0 / dg + 1.0)][4];
		
		// get the equation definition
		equation = getModelEquation (modelNumber);
		var = getModelVariables (modelNumber);
		x = getModelIndependentVar (modelNumber);
		
		// form the equation fx as ExpNode
		f = Matrix.infixToExpNode(equation, var);
		
		Bifurcation.setPitchfork(kValue, 0, "k", "g", 0.5, 2, 0.001, f, var,
				new double [] {0.0, 0.0, 1.25, 0.8, 0.0, 0.0}, x, "re", "im");
		Bifurcation.setPitchfork(kValue, 0, "k", "g", 0.5, 2, -0.001, f, var,
				new double [] {0.0, 0.0, 1.25, 0.8, 0.0, 0.0}, x, "re", "im");
		Bifurcation.setAndronov(kValue, 1, "k", "g", 0, 1, 0.001, f, var,
				new double [] {0.0, 0.0, 0.8, 0.8, 0.0, 0.6}, x, "re", "im");
		Bifurcation.setAndronov(kValue, 1, "k", "g", 0, 1, -0.001, f, var,
				new double [] {0.0, 0.0, 0.8, 0.8, 0.0, 0.6}, x, "re", "im");
		Bifurcation.setAndronov(kValue, 2, "k", "g", 1, 2, 0.001, f, var,
				new double [] {0.5546, 0.6932, 1.539, 0.8, 0.0, 0.599996}, 
				x, "re", "im");
		Bifurcation.setAndronov(kValue, 2, "k", "g", 1, 2, -0.001, f, var,
				new double [] {0.5546, 0.6932, 1.539, 0.8, 0.0, 0.599996}, 
				x, "re", "im");
		Bifurcation.setTangentBifLCycle(kValue, 3, "k", "g", 1.01, 2, -0.001, f, 
				new String [] {"x", "y", "g", "k", "re", "im"},
				new double [] {0.3, 0.0, 1.4, 0.87, 1.0, 0.0},
				//new double [] {0.5705838688419126, 0.0, 2.0, 0.77594, 1.0, 0.0},
				new String [] {"x", "y"}, 1, "re", "im");
		
		// print output
		setFile ();
		try {
			// set the header
			mWriter.write("g\tk1\tk2\tk3");
			mWriter.flush();
			
			// print result to file
			for (int i = 0; i < kValue.length; ++i) {
				mWriter.write(String.format ("\r\n%f", (i * dg)));
				for (int j = 0; j < kValue[i].length; ++j)
					mWriter.write(String.format("\t%f", kValue[i][j]));
				mWriter.flush();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	// method to return the model equation
	private StringBuffer [] getModelEquation (int modelName)
	{
		switch (modelName) {
		case 1:		// Model 17
			return new StringBuffer [] {
					new StringBuffer ("-y + tanh (g * x)"),
					new StringBuffer ("x - k * y")
			};
			
		case 2:		// Ueta2004
			return new StringBuffer [] {
					new StringBuffer ("-y1 + tanh (g * x1) - d * (x1 - x2)"),
					new StringBuffer ("x1 - k1 * y1"),
					new StringBuffer ("-y2 + tanh (g * x2) - d * (x2 - x1)"),
					new StringBuffer ("x2 - k2 * y2")
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
			return new String [] {"x", "y", "g", "k", "re", "im"};
		case 2:
			return new String [] {"x1", "y1", "x2", "y2", "g", "k1", "k2", "d",
					"re", "im"};
		case 3:
			return new String [] {"x1", "x2", "x3", "x4", "g", "k", "d", "re", "im"};
		case 4:
			return new String [] {"x1", "x2", "x3", "x4", "x5", "g1", "g2",
					"k1", "k2", "d", "s", "e1", "e3", "re", "im"};
		default:
			return null;
		}
	}
	
	// method to return the independent variables
	private String [] getModelIndependentVar (int modelName)
	{
		switch (modelName) {
		case 1: // model 17
			return new String [] {"x", "y"};
		case 2: case 3:
			return new String [] {"x1", "y1", "x2", "y2"};
		case 4:
			return new String [] {"x1", "x2", "x3", "x4", "x5"};
		default:
			return null;
		}
	}
	
	// method that sets the file writer
	private void setFile ()
	{
		// instantiate file and its writer
		try {
			int count = 0;
			java.io.File file = new java.io.File ("C:\\Users\\Adesola\\Desktop\\" +
					"m" + modelNumber + "contBif.txt");
			while (file.exists()) {
				file = new java.io.File (
						"C:\\Users\\Adesola\\Desktop\\m" + modelNumber +
						"contBif" + ++count + ".txt");
			}
			file.createNewFile();
			mWriter = new java.io.FileWriter(file);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new ContinuityTest ();
	}

}
