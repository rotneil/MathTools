package math.tools;

import java.io.IOException;

/**
 * This is an interface class for analyzing the bifurcation pattern of
 * a nonlinear dynamical system. It contains static methods for performing
 * each of these operations
 * 
 * @author Nehemiah Oluwafemi
 *
 */
public final class Bifurcation {
	
	/**
	 * This method analyzes a range of limit cycles along a straight line and prints
	 * the stability of the limit cycles for observation. The iteration process 
	 * starts from the value given in param0 to paramf
	 *  and loops until a bifurcation point is found or an equilibrium is 
	 * found.
	 * 
	 * @param file The file into which the stability values are printed
	 * @param f The model equation
	 * @param x The independent variables of the state equation
	 * @param variables The variables in the model
	 * @param values The corresponding value of the variables
	 * @param controlParameter This is the testing variable along which the test is taken
	 * @param param0 The starting point for the independent parameter
	 * @param paramf The end point for the independent parameter
	 * @param sectionIndex The Poincare section index in x
	 * @param dp The incremental value for the control parameter
	 * 
	 * @throws IllegalStateException If an equilibrium point is found.
	 */
	public static void showLimitCycleStability (java.io.File file, ExpNode [] f,
			String [] x, String [] variables, double [] values, 
			String controlParameter, double param0, double paramf, 
			int sectionIndex, double dp)
	{
		// initialize the number of the independent variables
		int n = x.length;
		double dt = 1E-2;
			
		// first define the time index and incorporate its variable name
		// and value
		int tIndex = variables.length;
		String t = "time";
		String [] var = new String [tIndex + 1];
		double [] val = new double [tIndex + 1];
		System.arraycopy(variables, 0, var, 0, variables.length);
		System.arraycopy(values, 0, val, 0, values.length);
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
		
		// set the file and the writer
		java.io.FileWriter writer = setFileWriter (null);
					
		// set the file details
		try {
			// set the title header
			String header = "Limit Cycle Stability test at ";
			String [] parameter = ExpNode.getParameters(x, variables);
			for (String str : parameter)
				header += str + " = " + ExpNode.getValue(str, variables, values)
				+ ", ";
			writer.write(header);
			
			// set the table header format
			header = "\r\n" + controlParameter + "\t\t";
			for (int i = 0; i < n; ++i)
				header += "h" + (i + 1) + "\t\t";
			for (String dx : Dx)
				header += dx + "\t\t";
			writer.write(header);
			
			writer.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		// start the analysis
		try {
			// loop through all the k-values
			for (double k = param0;	k >= param0 && k <= paramf; k += dp) {
				val[ExpNode.getIndex(controlParameter, var)] = k;
				
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
				
				// write the content of u to file
				writeToFile(writer, roots, var, val, DxIndices, controlParameter, dp);
			}
		} catch (IllegalStateException e) {
			try {
				writer.write("\r\nStable Equilibrium point encountered");
				writer.flush();
			} catch (IOException ee) {}
		}
	}
	
	
	/**
	 * This method analyzes a range of k-values for the bifurcation point
	 * of the limit cycle stability and prints to the give file all the 
	 * computed stability for this range. The iteration process starts from the value given in 
	 * values and loops until a bifurcation point is found or an equilibrium is 
	 * found.
	 * 
	 * @param f The model equation
	 * @param bifVar The bifurcation parameter
	 * @param bif0 The initial value of the bifurcation param
	 * @param bifMax The final testable value of the bifurcation param
	 * @param x The independent variables of the state equation
	 * @param variables The variables in the model
	 * @param values The corresponding value of the variables
	 * @param sectionIndex The Poincare section index in x
	 * @param dk k-value interval
	 * 
	 * @return the values of the first variational variables
	 * 
	 * @throws IllegalStateException If an equilibrium point is found.
	 */
	public static double [] getLimitCycleBifPoint (ExpNode [] f, String bifVar,
			double bifMax, String [] x, String [] variables, double [] values, 
			int sectionIndex, String re, String im, double dk)
	{
		// initialize the number of the independent variables
		int n = x.length;
		double dt = 1E-2;
		
		// first define the time index and incorporate its variable name
		// and value
		int tIndex = variables.length;
		String t = "time";
		String [] var = new String [tIndex + 1];
		double [] val = new double [tIndex + 1];
		System.arraycopy(variables, 0, var, 0, variables.length);
		System.arraycopy(values, 0, val, 0, values.length);
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
		
		// define the possible return array
		double [] mVal = new double [n + 1];
		
		// start the analysis
		try {
			// define the starting value for the bifurcation parameter
			double bif0 = ExpNode.getValue(bifVar, var, val);
			
			// loop through all the k-values
			for (double k = ExpNode.getValue(bifVar, var, val);
					k >= bif0 && k <= bifMax; k += dk) {
				
				// save the last computed values.
				// This is the most probable result.
				for (int i = 0; i < n; ++i)
					mVal [i] = val[DxIndices[i]];
				mVal[n] = ExpNode.getValue(bifVar, var, val1);
				
				// assign the new k to the val
				val[ExpNode.getIndex(bifVar, var)] = k;
				val = Algebra.getNextCrossing(f, xIndices, 
							sectionIndex, var, val, tIndex, dt, 1.0E-18);
				
				// reset the correction iteration error
				error = 1.0;
				while (error > 1E-9) {
					try {
						// initialize the variables of the variational equation
						val1 = Matrix.getVariationInitVal (val, x);
						
						// get the next poincare section crossing
						val1 = Algebra.getNextCrossing(fx, xIndices, 
								sectionIndex, var, val1, tIndex, dt, 1.0E-12);
						
						// set the value of xx
						for (int i = 0; i < DxIndices.length; ++i)
							xx[i] = val[DxIndices[i]];
						
						// reset ffx
						for (int i = 0; i < n - 1; ++i)
							ffx [i] = val1 [DxIndices[i]] - val[DxIndices[i]];
						ffx[n - 1] = val[xIndices[sectionIndex]];
						
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
			}
		} catch (IllegalStateException e) {}
		
		return mVal;
	}
	
	/**
	 * This method computes the tangent bifurcation pattern for a nonlinear system
	 * represented by f for the independent variables from x0 to xf at an interval 
	 * dx. The bifurcation values are put in the appropriate indices of array
	 * bif through the indicated index. The info, if specified is the additional
	 * value s
	 * 
	 * @param bif The array for all the bifurcation patterns
	 * @param index The column in bif that this bifurcation patterns should be
	 * saved.
	 * @param kAxis The value on the vertical axis of the 2-param bifurcation diagram
	 * @param gAxis The horizontal axis for the 2-param bifurcation diagram
	 * @param g0 The starting point of the independent param
	 * @param gf The ending point of the independent param
	 * @param dg The interval of evaluation
	 * @param f The state equation of the nonlinear system
	 * @param variables The variables defined in the ExpNode f which should include 
	 * the variables for the real and imaginary part of the stability eigenvalue.
	 * @param values The initial values of the variables which should include 
	 * the values for the real and imaginary part of the stability eigenvalue
	 * @param x The independent variables of f
	 * @param sectionIndex The index of the independent coordinate 
	 * (in x above) about which the Poicare section is taken
	 * @param re The real part of the complex eigenvalue usually zero for 
	 * pitchfork bifurcation of the equilibrium.
	 * @param im The imaginary part of the complex eigenvalue usually zero for a 
	 * pitchfork bifurcation of the equilibrium
	 */
	public static void setTangentBifLCycle (double [] [] bif, int index, 
			String kAxis, String gAxis, double g0, double gf, double dg, 
			ExpNode [] f, String [] variables, double [] values, String [] x, 
			int sectionIndex, String re, String im)
	{	
		// initialize the number of the independent variables
		int n = x.length;
		double dt = 1E-2;
		
		// first define the time index and incorporate its variable name
		// and value
		int tIndex = variables.length;
		String t = "time";
		String [] var = new String [tIndex + 1];
		double [] val = new double [tIndex + 1];
		System.arraycopy(variables, 0, var, 0, tIndex);
		System.arraycopy(values, 0, val, 0, tIndex);
		var [tIndex] = t;
		val [tIndex] = 0.0;
		
		// redefine variables to suit second variation equation
		var = Matrix.getSecondVariationVar(var, x);
		double [] val1 = Matrix.getSecondVariationInitVal (val, x);	// un
		ExpNode [] f2x = Matrix.getSecondVariationEq(f, x, kAxis);
		
		// put the embedded variables into fx
		int [] [] Dhx = Matrix.Dhx(x, sectionIndex);
		int [] [] Dhu = Matrix.Dhu(x, sectionIndex);
		int [] Dqx = Matrix.Dqx (x, sectionIndex);
		
		// form the derivative function for the Newton-Raphson iteration
		ExpNode [] [] dw = Matrix.getFirstVariationXNode(x);
		ExpNode [] dk = Matrix.getFirstVariationParamNode(x);
		ExpNode [] [] DFx = new ExpNode [n + 1] [n + 1];
		
		// first row element
		ExpNode [] [] DTu = Matrix.product(Dhx, Matrix.product(dw, Dhu));
		ExpNode [] DTt = Matrix.product(Dhx, f);
		ExpNode [] DTk = Matrix.product(Dhx, dk);
		
		// complete DFx for the first row
		for (int i = 0; i < n - 1; ++i) {
			for (int j = 0; j < n - 1; ++j)
				// subtract the identity element from diagonal elements of DTu
				DFx [i][j] = (i == j ? BinOpNode.make ('-', DTu [i][j], 
						ExpNode.makeConstantNode(1.0)) : DTu [i][j]);
			DFx [i][n - 1] = DTt[i];
			DFx [i][n] = DTk [i];
		}
		
		// second row element
		ExpNode [] Dqu = Matrix.product(Dqx, Matrix.product(dw, Dhu));
		ExpNode Dqt = Matrix.product(Dqx, f);
		ExpNode Dqk = Matrix.product(Dqx, dk);
		
		// complete second row
		for (int i = 0; i < n - 1; ++i)
			DFx [n - 1][i] = Dqu[i];
		DFx [n - 1] [n - 1] = Dqt;
		DFx [n - 1] [n] = Dqk;
		
		// third row elements and assign the derivative of X wrt to u as the value of X
		ExpNode X = Matrix.getX(x, re, im, kAxis);
		ExpNode duX = X.derivative("re");
		ExpNode [] dXx = Matrix.derivate(X, x);
		ExpNode [] DXu = Matrix.product(dXx, Dhu);
		ExpNode DXt = Matrix.product(dXx, f);
		ExpNode DXk = X.derivative(kAxis);
		//ExpNode DXk = Matrix.product(dXx, dk);
		
		// complete third row after differentiating wrt u=re
		for (int i = 0; i < n - 1; ++i)
			DFx[n][i] = DXu[i].derivative("re");
		DFx[n][n - 1] = DXt.derivative("re");
		DFx[n][n] = DXk.derivative("re");
		
		
		// Newton-raphson iteration
		int count = 0;
		double error = 1.0;
		String [] h = Matrix.hx(x, sectionIndex);
		String [] Dx = new String [n + 1];
		int ttIndex = n - 1;
		for (int i = 0; i < h.length; ++i) 
			Dx[i] = h[i];
		Dx[ttIndex] = t;
		Dx[n] = kAxis;
		
		int [] DxIndices = ExpNode.getXIndex(Dx, var);
		int [] xIndices = ExpNode.getXIndex(
				Matrix.getSecondVariationIndependentVar(x), var);
		double [] xx = new double [n + 1];
		for (int i = 0; i < xx.length; ++i)
			xx [i] = val[DxIndices[i]];
		double [] ffx = new double [n + 1];
		double [] [] Dff = new double [n + 1][n + 1];
		
		// define and initialize the last stable values
		double [] stableVal = new double [n + 1];
		for (int i = 0; i < n + 1; ++i)
			stableVal [i] = val [DxIndices[i]];
		
		// begin the iteration
		for (double gg = ExpNode.getValue(gAxis, var, val); 
				gg >= g0 && gg <= gf; gg += dg) {
			// set the iteration value of the independent coord
			val[ExpNode.getIndex(gAxis, variables)] = gg;
			val = Algebra.getNextCrossing(f, xIndices, 
					sectionIndex, var, val, tIndex, dt, 1.0E-18);
			
			// reset the correction iteration error
			error = 1.0;
			count = 0;
			while (error > 1E-6) {
				try {
					// initialize the variables of the variational equation
					val1 = Matrix.getSecondVariationInitVal (val, x);
					
					// get the next poincare section crossing
					val1 = Algebra.getNextCrossing(f2x, xIndices, 
							sectionIndex, var, val1, tIndex, dt, 1.0E-18);
					
					// initialize xx
					for (int i = 0; i < DxIndices.length; ++i)
						xx[i] = val[DxIndices[i]];
					
					// reset ffx
					for (int i = 0; i < n - 1; ++i)
						ffx [i] = val1 [DxIndices[i]] - val[DxIndices[i]];
					ffx[n - 1] = val[xIndices[sectionIndex]];
					ffx[n] = duX.value(var, val1);
					
					Dff = Matrix.values(DFx, var, val1);
					
					ffx = Matrix.LUSolver(Dff, 
							Matrix.subtract(Matrix.product(Dff, xx), ffx));
					error = Matrix.getMaxError(ffx, xx);
					
					// reset the elements of Dx in val
					for (int i = 0; i < Dx.length; ++i)
						val[DxIndices[i]] = ffx[i];
					
					// check for a certain number of iterations
					if (++count > 200) {
						// adjust the val to the last stable values
						for (int i = 0; i < n + 1; ++i)
							val[DxIndices[i]] = stableVal[i];
						
						val1 = Bifurcation.getLimitCycleBifPoint(f, kAxis, 1.0, x, 
								variables, val, sectionIndex, re, im, 0.00001);
						// reset the elements of Dx in val
						for (int i = 0; i < Dx.length; ++i)
							val[DxIndices[i]] = val1[i];
						
						break;
					}
				} catch (MatrixException e) {
					break;
				} catch (IllegalStateException e) {
					// at this state, the orbit has gone into the equilibrium
					// position, choose a convenient point outside the limit cycle
					// adjust the val to the last stable values
					for (int i = 0; i < n + 1; ++i)
						val[DxIndices[i]] = stableVal[i];
					
					val1 = Bifurcation.getLimitCycleBifPoint(f, kAxis, 1.0, x, 
							variables, val, sectionIndex, re, im, 0.00001);
					// reset the elements of Dx in val
					for (int i = 0; i < Dx.length; ++i)
						val[DxIndices[i]] = val1[i];
					
					break;
				}
			}
			
			System.out.printf ("%f\t%f\t%f\n", gg, val[DxIndices[n]], 
					val[DxIndices[0]]);
			
			// save the convergence values into the array bif
			for (int i = 0; i < n + 1; ++i)
				stableVal [i] = val[DxIndices[i]];
			
			bif [getIndex(gg, dg)] [index] = ExpNode.getValue(kAxis, var, val);
		}
	}
	
	/**
	 * This method computes pitchfork bifurcation pattern for a nonlinear system
	 * represented by f for the independent variables from x0 to xf at an interval 
	 * dx. The bifurcation values are put in the appropriate indices of array
	 * bif through the indicated index. The info, if specified is the additional
	 * value s
	 * 
	 * @param bif The array for all the bifurcation patterns
	 * @param index The column in bif that this bifurcation patterns should be
	 * saved.
	 * @param kAxis The value on the vertical axis of the 2-param bifurcation diagram
	 * @param gAxis The horizontal axis for the 2-param bifurcation diagram
	 * @param g0 The starting point of the independent param
	 * @param gf The ending point of the independent param
	 * @param dg The interval of evaluation
	 * @param f The state equation of the nonlinear system
	 * @param var The variables defined in the ExpNode f
	 * @param values The initial values of the variables
	 * @param x The independent variables of f
	 * @param re The real part of the complex eigenvalue usually zero for 
	 * pitchfork bifurcation of the equilibrium.
	 * @param im The imaginary part of the complex eigenvalue usually zero for a 
	 * pitchfork bifurcation of the equilibrium
	 */
	public static void setPitchfork (double [] [] bif, int index, String kAxis, 
			String gAxis, double g0, double gf, double dg, ExpNode [] f, 
			String [] var, double [] values, String [] x, String re, String im)
	{
		// redefine another values
		double [] val = new double [values.length];
		for (int i = 0; i < values.length; ++i)
			val [i] = values[i];
		
		// from f form the stability function
		ExpNode X = Matrix.getX(f, x, re, im);
		
		// for the equation of Newton-Raphson interpol
		ExpNode [] fx = new ExpNode [f.length + 1];
		for (int i = 0; i < f.length; ++i)
			fx [i] = f[i];
		fx[f.length] = X;
		
		// for the independent variable for the newton-raphson such that
		// it includes the continuity variable "k"
		String [] Dx = new String [x.length + 1];
		for (int i = 0; i < x.length; ++i)
			Dx[i] = x[i];
		Dx[x.length] = kAxis;
		
		// differentiate wrt to Dx
		ExpNode [] [] Dfx = Matrix.jacobian(fx, Dx);
		
		// perform either forward or backward analysis
		double error = 1.0;
		int [] DxIndices = ExpNode.getXIndex(Dx, var);
		double [] xx = new double [DxIndices.length];
		for (int i = 0; i < xx.length; ++i)
			xx[i] = val[DxIndices[i]];
		double [] ffx = new double [fx.length];
		double [] [] Dff = new double [xx.length][xx.length];
		
		// begin the iteration
		for (double gg = ExpNode.getValue(gAxis, var, val); 
				gg >= g0 && gg <= gf; gg += dg) {
			// set the iteration value of the independent coord
			val[ExpNode.getIndex(gAxis, var)] = gg;
			
			// reset the correction iteration error
			error = 1.0;
			while (error > 1E-6) {
				try {
					ffx = Matrix.values(fx, var, val);
					Dff = Matrix.values(Dfx, var, val);
					ffx = Matrix.LUSolver(Dff, Matrix.subtract(
							Matrix.product(Dff, xx), ffx));
					error = Matrix.getMaxError(ffx, xx);
					xx = ffx;
					
					// return the new values into val
					for (int i = 0; i < xx.length; ++i)
						val[DxIndices[i]] = xx[i];
				} catch (MatrixException e) {
					break;
				} 
			}
			
			// save the convergence values into the array bif
			bif [getIndex(gg, dg)] [index] = val[ExpNode.getIndex(kAxis, var)];
		}
	}
	
	/**
	 * This method computes the bifurcation parameter values for andronov 
	 * bifurcation of the equilibrium point. The values are put at the 
	 * appropriate index in 2-dimensional array bif using int index and the
	 * calculation of the index in bif[index].
	 * @param bif The 2-dimensional array that keeps the bifurcation values
	 * @param index The index in bif where the bifurcation values are kept
	 * @param kAxis The String name for the vertical axis of the 2-parameter
	 * bifurcation plot
	 * @param gAxis The horizontal axis name for the 2-parameter bifurcation plot
	 * @param g0 The starting point on the horizontal axis
	 * @param gf The end point on the horizontal axis
	 * @param dg The interval for the horizontal axis
	 * @param f The state functions as an array of ExpNode
	 * @param var The variables of the state functions
	 * @param values The corresponding values of the variables in var
	 * @param x The independent coordinates of the state function
	 * @param re The real part of the stability eigenvalue
	 * @param im The imaginary part of the complex stability eigenvalue
	 */
	public static void setAndronov (double [] [] bif, int index, String kAxis, 
			String gAxis, double g0, double gf, double dg, ExpNode [] f, 
			String [] var, double [] values, String [] x, String re, String im)
	{
		// redefine another values
		double [] val = new double [values.length];
		for (int i = 0; i < values.length; ++i)
			val [i] = values[i];
		
		// from f form the stability function
		ExpNode X = Matrix.getX(f, x, re, im);
		
		// for the equation of Newton-Raphson interpol
		ExpNode [] fx = new ExpNode [f.length + 2];
		for (int i = 0; i < f.length; ++i)
			fx [i] = f[i];
		fx [f.length] = X;
		fx [f.length + 1] = X;
		
		// for the independent variable for the newton-raphson such that
		// it includes the continuity variable "k"
		String [] Dx = new String [x.length + 2];
		for (int i = 0; i < x.length; ++i)
			Dx[i] = x[i];
		Dx[x.length] = kAxis;
		Dx[x.length + 1] = im;
		
		// differentiate wrt to Dx
		ExpNode [] [] Dfx = Matrix.jacobian(fx, Dx);
		
		// perform either forward or backward analysis
		double error = 1.0;
		int [] DxIndices = ExpNode.getXIndex(Dx, var);
		double [] xx = new double [DxIndices.length];
		for (int i = 0; i < xx.length; ++i)
			xx[i] = val[DxIndices[i]];
		double [] ffx = new double [fx.length];
		double [] [] Dff = new double [fx.length][fx.length];
		
		// begin the iteration
		for (double gg = ExpNode.getValue(gAxis, var, val); 
				gg >= g0 && gg <= gf; gg += dg) {
			// set the iteration value of the independent coord
			val[ExpNode.getIndex(gAxis, var)] = gg;
			
			// reset the correction iteration error
			error = 1.0;
			while (error > 1E-12) {
				try {
					ffx = Matrix.values(fx, var, val, x.length, x.length + 1);
					Dff = Matrix.values(Dfx, var, val, x.length, x.length + 1);
					ffx = Matrix.LUSolver(Dff, Matrix.subtract(
							Matrix.product(Dff, xx), ffx));
					error = Matrix.getMaxError(ffx, xx);
					xx = ffx;
					
					// return the new values into val
					for (int i = 0; i < xx.length; ++i)
						val[DxIndices[i]] = xx[i];
				} catch (MatrixException e) {
					break;
				} 
			}
			
			// save the convergence values into the array bif
			bif [getIndex(gg, dg)] [index] = val[ExpNode.getIndex(kAxis, var)];
		}
	}
	
	private static int getIndex (double gg, double dg) { 
		return ((int) (gg / Math.abs(dg) + 0.5));
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
	 * @param DxIndices The indices of the Poincare variables
	 * @param var The variables of the state model
	 * @param val The corresponding values of the state variables
	 * @param dk The incremental vertical value
	 */
	private static void writeToFile (java.io.FileWriter writer, 
			ComplexNumber [] eigen,	String [] var, double [] val, 
			int [] DxIndices, String controlParameter, double dk) 
	{
		// begin writing into the file
		try {
			// print result to file
			writer.write(String.format ("\r\n%f", 
					val[ExpNode.getIndex(controlParameter, var)]));
			for (int j = 0; j < eigen.length; ++j)
				writer.write(String.format("\t%s", eigen[j].format()));
			for (int j = 0; j < DxIndices.length; ++j)
				writer.write(String.format("\t%f", val[DxIndices[j]]));
			
			writer.flush();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
}
