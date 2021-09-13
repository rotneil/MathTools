package phDModels;

import java.awt.Color;
import java.awt.FileDialog;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.text.DecimalFormat;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.JOptionPane;

import math.tools.Algebra;
import math.tools.ExpNode;
import math.tools.Matrix;

public class PhasePlotter extends javax.swing.JFrame
{
	private static final long serialVersionUID = 1L;
	
	// instance variables
	private static int n = 100000;
	private double dt = 0.01;
	
	private double x0 = -3.0;
	private double xf = 3.0;
	private double y0 = -3.0;
	private double yf = 3.0;
	
	private static int xStart = 50;
	private static int yStart = 50;
	private static int width = 400;
	private static int height = 400;
	
	// State model variables
	private java.util.ArrayList<PlotModel> plotModel;
	private int currentModel;
	private String modelName;
	private String fileName = new String ("modelFile.txt");
	private ExpNode [] fx;
	private String [] x;
	private int [] xIndices;
	private int [] parIndices;
	private int [] controlParam;
	private double dp = 0.0001;
	private String [] var;
	private String [] par;
	private double [] val;
	private String [] plotVar;
	private int [] plotIndices;
	private double [] [] xVal;
	private double [] sx, sy;
	
	private boolean multiplePlot = false;
	private boolean showSteadyState = false;
	private int steadyStateCounter = 1;
	private int plotCounter = n;
	private boolean iteratedPlot = false;
	private boolean backward = false;
	
	private boolean specialPlot = false;
	private String specialX = "", specialY = "";
	private String specialLegendX = "", specialLegendY = "";
	
	// method to initialize the important plot variables
	private void initializePlotter ()
	{
		// initialize the file
		java.io.File modelFile = new java.io.File (fileName);
		
		// if the file exists, upload plot variables
		if (modelFile.exists())
			uploadModelFromFile (modelFile);
		else {
			modelName = "BVP Oscillator";
			setTitle (modelName);
			x = new String [] {"x", "y"};
			var = new String [] {"g", "k", "x", "y"};
			val = new double [] {1.5, 0.8, 2.5, 0.0};
			fx = Matrix.infixToExpNode(new StringBuffer [] {
					new StringBuffer ("-y + tanh (g * x)"),
					new StringBuffer ("x - k * y")}, var);
			xIndices = ExpNode.getXIndex(x, var);
			par = ExpNode.getParameters(x, var);
			parIndices = ExpNode.getXIndex(par, var);
			plotVar = new String [] {"x", "y"};
			plotIndices = ExpNode.getXIndex(plotVar, x);
			controlParam = new int [] {parIndices[0], parIndices[1]};
			xVal = new double [n][x.length];
			sx = new double [n];
			sy = new double [n];
			plotCounter = n;
			
			// instantiate plot models
			plotModel = new java.util.ArrayList<>();
			plotModel.add(new PlotModel (modelName, fx, x, var, val, plotIndices));
			currentModel = 0;
		}
	}
	
    // method to upload a model from file
	private void uploadModelFromFile (java.io.File modelFile)
	{
		try {
			// local variables
			String modelName;
			ExpNode [] fx;
			String [] eq;
			String [] x;
			String [] var;
			String [] valStr;
			double [] val;
			String [] cStr;
			int [] coord;
			
			// open file, if it exists
			java.util.Scanner scanner = new java.util.Scanner(modelFile);
			scanner.useDelimiter("\r\n");
			
			// instantiate plot models to zero
			plotModel = new java.util.ArrayList<>();
			
			// read current model from file
			currentModel = Integer.parseInt(scanner.nextLine());
			
			// start reading the content of the file
			while (!scanner.nextLine().equals("END")) {
				// read the model name
				modelName = scanner.nextLine();
				
				// read the equation, variables and x before converting to ExpNode
				eq = scanner.nextLine().split(",");
				x = scanner.nextLine().split(",");
				var = scanner.nextLine().split(",");
				valStr = scanner.nextLine().split(",");
				cStr = scanner.nextLine().split(",");
				fx = new ExpNode [eq.length];
				for (int i = 0; i < eq.length; ++i)
					fx [i] = Algebra.infixToBinOpNode(
							new StringBuffer (eq[i]), var);
				
				val = new double [valStr.length];
				for (int i = 0; i < val.length; ++i)
					val[i] = Double.parseDouble(valStr[i]);
				
				coord = new int [cStr.length];
				for (int i = 0; i < coord.length; ++i)
					coord[i] = Integer.parseInt(cStr[i]);
				
				// include the read model into plot models
				plotModel.add(new PlotModel (modelName, fx, x, var, val, coord));
			}
			setCurrentModel ();
			scanner.close();
		} catch (NumberFormatException e) {
			javax.swing.JOptionPane.showMessageDialog(this, e.getMessage(), 
					e.getClass().getSimpleName(), 
					javax.swing.JOptionPane.ERROR_MESSAGE);
		} catch (IllegalStateException e) {
			javax.swing.JOptionPane.showMessageDialog(this, e.getMessage(), 
					e.getClass().getSimpleName(), 
					javax.swing.JOptionPane.ERROR_MESSAGE);
		} catch (java.io.FileNotFoundException e) {
			javax.swing.JOptionPane.showMessageDialog(this, e.getMessage(), 
					e.getClass().getSimpleName(), 
					javax.swing.JOptionPane.ERROR_MESSAGE);
		}
	}
	
	// method to allow user to select model
	private void selectModel () 
	{
		try {
			// prompt user to select for the list of models
			PlotModel model = (PlotModel) javax.swing.JOptionPane.showInputDialog(
					this, "Select model from list: ", "Model Selection Dialog",
					javax.swing.JOptionPane.INFORMATION_MESSAGE,
					null, plotModel.toArray(), plotModel.get(0));
			
			// get model's index
			for (int i = 0; i < plotModel.size(); ++i) {
				if (model.equals(plotModel.get(i))) {
					currentModel = i;
					setCurrentModel ();
					getSolution ();
					writePlotModelToFile();
					return;
				}
			}
		} catch (Exception e) {
			javax.swing.JOptionPane.showMessageDialog(this, e.getMessage(), 
					e.getClass().getSimpleName(), 
					javax.swing.JOptionPane.ERROR_MESSAGE);
		}
	}
	
	// method to set the values of the parameter
	private void setParameter  ()
	{
		try {
			String prompt = "Set the values for parameters as: \n";
			for (String param : par)
				prompt += param + ", ";
			prompt = prompt.substring(0, prompt.length() - 2);
			
			String input = javax.swing.JOptionPane.showInputDialog(this,
					prompt, "Parameter Value Dialog",
					javax.swing.JOptionPane.INFORMATION_MESSAGE);
			
			if (input == null)
				throw new NullPointerException ();
			
			// split input into values
			String [] value = input.split(",");
			
			// check for out of bounds exception
			if (value.length != parIndices.length)
				throw new IllegalArgumentException ("Insufficient entry! " +
						"Check that each parameter is given a value.");
			
			// set the double values
			for (int i = 0; i < parIndices.length; ++i)
				val[parIndices[i]] = Double.parseDouble(value[i]);
			
			getSolution ();
			showParameters();
			writePlotModelToFile ();
		} catch (NumberFormatException e) {
			javax.swing.JOptionPane.showMessageDialog(this, e.getMessage(), 
					e.getClass().getSimpleName(), 
					javax.swing.JOptionPane.ERROR_MESSAGE);
		} catch (NullPointerException e) {}
	}
	
	// method to change the plot parameters
	private void changeParameter ()
	{
		try {
			String prompt = "Choose the two comma separated " +
					"control parameters from list: ";
			for (String p : par)
				prompt += p + ", ";
			prompt = prompt.substring(0, prompt.length() - 2);
			
			String input = javax.swing.JOptionPane.showInputDialog(this, 
					prompt, "Control Parameter Dialog", 
					javax.swing.JOptionPane.INFORMATION_MESSAGE);
			String [] parStr = input.split(",");
			
			// check and confirm valid input
			if (parStr.length != 2)
				throw new IllegalArgumentException ("The number of entered " +
						"parameters can only be two");
			
			controlParam = ExpNode.getXIndex(parStr, var);
		} catch (IllegalArgumentException e) {
			javax.swing.JOptionPane.showMessageDialog(this, e.getMessage(), 
					e.getClass().getSimpleName(), 
					javax.swing.JOptionPane.ERROR_MESSAGE);
		} catch (NullPointerException e) {}
	}
	
	// method to set the current model
	private void setCurrentModel ()
	{
		PlotModel model = plotModel.get(currentModel);
		setTitle (model.modelName);
		x = model.x;
		var = model.var;
		val = model.val;
		fx = model.fx;
		xIndices = ExpNode.getXIndex(x, var);
		par = ExpNode.getParameters(x, var);
		parIndices = ExpNode.getXIndex(par, var);
		controlParam = new int [] {parIndices[0], parIndices[1]};
		plotIndices = model.coord;
		plotVar = new String [plotIndices.length];
		for (int i = 0; i < plotVar.length; ++i)
			plotVar[i] = x[i];
		xVal = new double [n][x.length];
		sx = new double [n];
		sy = new double [n];
		
		// set the plot initial values
		for (int i = 0; i < xIndices.length; ++i) 
			xVal[0][i] = val[xIndices[i]];
		sx[0] = xVal[0][plotIndices[0]];
		sy[0] = xVal[0][plotIndices[1]];
	}
	
	// method to show the model dialog and to initialize other variables
	private void showModelDialog ()
	{
		// instantiate the model dialog and use inner class annotation
		ModelDialog dialog = new ModelDialog (this, 
				new ModelDialog.ModelDialogListener()
		{
			@Override
			public void onModelInitialized(String model, ExpNode[] f, 
					String[] xStr, String[] variables, double[] values, 
					int [] coord) 
			{
				// include the current model to list and write to file
				plotModel.add(new PlotModel (model, f, xStr, variables,
						values, coord));
				currentModel = plotModel.size() - 1;
				writePlotModelToFile ();
				getSolution ();
				
			}
		});
		dialog.setVisible(true);
	}
	
	// method to write plot models to file
	private void writePlotModelToFile ()
	{
		// instantiate a third-party thread to write to file
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {
			public void run () {
				// TODO Auto-generated method stub
				try {
					// create model file if does not exit
					java.io.File modelFile = new java.io.File (fileName);
					modelFile.createNewFile();
					
					// instantiate the file writer
					java.io.FileWriter writer = new java.io.FileWriter (
							modelFile);
					
					// write the current mode to file
					writer.write("" + currentModel + "\r\n");
					
					// write each model to file
					for (PlotModel current : plotModel) {
						// write the model name
						writer.write("\r\n" + current.modelName);
						
						// write the model equation
						String eq = current.fx [0].toString();
						for (int i = 1; i < current.fx.length; ++i)
							eq += "," + current.fx[i].toString();
						writer.write("\r\n" + eq);
						
						// write the independent variables
						String xx = current.x [0];
						for (int i = 1; i < current.x.length; ++i)
							xx += "," + current.x[i];
						writer.write("\r\n" + xx);
						
						// write the Variables
						String var = current.var[0];
						for (int i = 1; i < current.var.length; ++i)
							var += "," + current.var[i];
						writer.write("\r\n" + var);
						
						// write the initial values
						String val = String.valueOf(current.val[0]);
						for (int i = 1; i < current.val.length; ++i)
							val += "," + current.val[i];
						writer.write("\r\n" + val);
						
						// write the plot coordinate
						String coord = String.valueOf(current.coord[0]);
						for (int i = 1; i < current.coord.length; ++i)
							coord += "," + current.coord[i];
						writer.write("\r\n" + coord);
						
						// separate models with empty line
						writer.write("\r\n");
					}
					
					// end file and close writer
					writer.write("END");
					writer.flush();
					writer.close();
				} catch (java.io.IOException e) {
					javax.swing.JOptionPane.showMessageDialog(PhasePlotter.this, 
							e.getMessage(), 
							e.getClass().getSimpleName(), 
							javax.swing.JOptionPane.ERROR_MESSAGE);
				}
			}
		});
		executor.shutdown();
	}
	
	// method to set the variable initial condition
	private void setPlotVariableInitial (double x, double y)
	{
		val [xIndices[plotIndices[0]]] = sx[0] = x;
		val [xIndices[plotIndices[1]]] = sy[0] = y;
		
		// set the initial condition
		for (int i = 0; i < xIndices.length; ++i)
			xVal[0][i] = val[xIndices[i]];
		
		getSolution ();
	}	// end method setVariableInitialCondition

	// method to toggle the iterated menu option
	private void setIteratedPlot ()
	{
		// toggle flag
		iteratedPlot = !iteratedPlot;
		
		// modify plot state
		if (iteratedPlot)
			plotCounter = steadyStateCounter;
		else
			plotCounter = xVal.length;
		
		// modify the menu
		iteratedMenu.setState(iteratedPlot);
		phasePlot.repaint();
	}
	
	// method to set the parameter increment value
	private void setParameterIncrement ()
	{
		try {
			dp = Double.parseDouble(
					javax.swing.JOptionPane.showInputDialog(this,
						"Enter the incremental value for control parameter:",
						"Input Dialog",
						javax.swing.JOptionPane.INFORMATION_MESSAGE));
		} catch (IllegalArgumentException e) {
			javax.swing.JOptionPane.showMessageDialog(this, 
					e.getClass().getSimpleName(), e.getMessage(),
					javax.swing.JOptionPane.ERROR_MESSAGE);
		} catch (NullPointerException e) {}
	}
	
	// method to set the steadyState iterate
	private void setSteadyState ()
	{
		// toggle showSteadyState
		showSteadyState = !showSteadyState;
		
		// prompt the user for entry
		if (showSteadyState) {
			try {
				String input = JOptionPane.showInputDialog(this, 
							"Enter the steady state iterate (less than " +
							n + ": ");
				// check for null input
				if (input == null)
					return;
				
				// convert the input to an integer
				steadyStateCounter = Integer.parseInt(input);
						
				// check for valid input for the steady state
				if (steadyStateCounter > n || steadyStateCounter < 1)
					throw new IllegalArgumentException ("Invalid entry for " +
							"the steady state iterate");
				
				steadyMenu.setState(showSteadyState);
				getSolution ();
			} catch (NumberFormatException nfe) {
				JOptionPane.showMessageDialog(this, nfe.getMessage(),
						nfe.getClass().getSimpleName(), JOptionPane.ERROR_MESSAGE);
				showSteadyState = false;
				steadyStateCounter = 1;
			} catch (NullPointerException npe) {}
		} else {
			steadyStateCounter = 1;
			steadyMenu.setState(showSteadyState);
			getSolution ();
		}
	}
	
    // method to set the initial condition for all the variables
	private void setVariablesInitial ()
	{
		// prompt user for input
		try {
			String variable = "";
			for (String v : x)
				variable += v + ", ";
			variable = variable.substring(0, variable.length() - 2);
			
			String input = javax.swing.JOptionPane.showInputDialog(this,
					"Enter a comma separated initial value for " + variable, 
					"Variable Initial Value",
					javax.swing.JOptionPane.INFORMATION_MESSAGE);
			String [] initial = input.split(",");
			
			// check for correct input
			if (initial.length != x.length)
				throw new IllegalArgumentException ("Invalid input. Check that " +
						"a value is given for each independent variable");
			
			// set the initial condition
			for (int i = 0; i < initial.length; ++i)
				val [xIndices[0]] = xVal [0][i] = Double.parseDouble(initial[i]);
			sx[0] = xVal [0][plotIndices[0]];
			sy[0] = xVal [0][plotIndices[1]];
			getSolution();
		} catch (IllegalArgumentException e) {
			javax.swing.JOptionPane.showMessageDialog(this, 
					e.getClass().getSimpleName(), e.getMessage(), 
					javax.swing.JOptionPane.INFORMATION_MESSAGE);
		} catch (NullPointerException e) {}
	}
	
    // method to solve for the system and plot 
    private void getSolution ()
    {
    	// iterate n-times
    	for (int i = 1; i < n; ++i) {
    		val = Algebra.rungeKutta(fx, xIndices, var, val, dt);
    		for (int j = 0; j < x.length; ++j)
    			xVal [i] [j] = val[xIndices[j]];
    		sx [i] = val[xIndices[plotIndices[0]]];
    		sy [i] = val[xIndices[plotIndices[1]]];
    	}
    	
    	// assign special plot variables
    	if (specialPlot)
    		setSpecialVariables ();
    	
    	// show the parameters
    	showParameters ();
    	phasePlot.repaint();
    }
    
    // method to assign variables for the special plot
    private void setSpecialVariables ()
    {
    	try {
    		// first assign for x
        	if (specialX.contains("+")) {
        		int index1 = ExpNode.getIndex(
        				specialX.substring(0, specialX.indexOf("+")), var);
        		int index2 = ExpNode.getIndex(specialX.substring(
        				specialX.indexOf("+"), specialX.length()), var);
        		int indexY = ExpNode.getIndex(specialY, var);
        		
        		// check the entered indices
        		int max = var.length;
        		if (index1 >= max || index2 >= max || indexY >= max)
        			throw new IllegalArgumentException (
        					"The entered index is out of range");
        		
        		// set the values for the special plot
        		for (int i = 0; i < n; ++i) {
        			sx[i] = xVal[i][index1] + xVal[i][index2];
        			sy[i] = xVal[i][indexY];
        		}
        		// set the legend for the special axis
        		specialLegendX = var[index1] + " + " + var[index2];
        		specialLegendY = var[indexY];
        	}
        	
        	// assign for y
        	if (specialY.contains("+")) {
        		int index1 = Integer.parseInt(
        				specialY.substring(0, specialY.indexOf("+")));
        		int index2 = Integer.parseInt(specialY.substring(
        				specialY.indexOf("+"), specialY.length()));
        		int indexX = Integer.parseInt(specialX);
        		
        		// check the entered indices
        		int max = var.length;
        		if (index1 >= max || index2 >= max || indexX >= max)
        			throw new IllegalArgumentException (
        					"The entered index is out of range");
        		
        		
        		// set the values for the special plot
        		for (int i = 0; i < sy.length; ++i) {
        			sy[i] = xVal[i][index1] + xVal[i][index2];
        			sx[i] = xVal[i][indexX];
        		}
        		
        		// set the legend for the special axis
        		specialLegendY = var[index1] + " + " + var[index2];
        		specialLegendX = var[indexX];
        	}
        	
    	} catch (IllegalArgumentException e) {
			javax.swing.JOptionPane.showMessageDialog(phasePlot, e.getMessage(),
					"INPUT ERROR", javax.swing.JOptionPane.ERROR_MESSAGE);
		}
    }
 	
	// method to change the plot variables
	private void changePlotVariables ()
	{
		// set the special plot to false
		specialPlot = false;
		
		try {
			// prepare the output message
			String varStr = x[0];
			for (int i = 0; i < x.length - 1; ++i)
				varStr += ", " + x[i];
			
			// obtain user's input
			String input = javax.swing.JOptionPane.showInputDialog(phasePlot, 
					"Enter the new variables from the list " + varStr
					+ " as comma separated input",
					javax.swing.JOptionPane.YES_NO_OPTION);
			
			// split input around comma
			String [] str = input.split(",");
			
			// confirm appropriate entry
			if (str.length != 2)
				throw new IllegalArgumentException ("Invalid input format!!!");
			
			// check for special plot
			if (str[0].contains("+") || str[1].contains("+")) {
				// flag special plot
				specialPlot = true;
				specialY = str[0];
				specialX = str[1];
			}
			else {
				// convert input
				int index1 = ExpNode.getIndex(str[0], x);
				int index2 = ExpNode.getIndex(str[1], x);
				
				// assign values
				plotIndices [0] = index1;
				plotIndices [1] = index2;
			}
			
			phasePlot.repaint();
		} catch (IllegalArgumentException e) {
			javax.swing.JOptionPane.showMessageDialog(phasePlot, e.getMessage(),
					"INPUT ERROR", javax.swing.JOptionPane.ERROR_MESSAGE);
		} catch (Exception e) {}
	}
	
	// method to change variable ranges
	private void changeVariableRange()
	{
		try {
			String input = JOptionPane.showInputDialog(phasePlot, "Enter the bounds as x0,xf,y0,yf",
					"Variable Range Alteration Dialog", JOptionPane.INFORMATION_MESSAGE);
			
			// split string around commas into parts
			String [] array = input.split(",");
			
			// confirm proper entry
			if (array.length != 4)
				throw new IllegalArgumentException ("The entered input does not " +
						"confirm to x0,xf,y0,yf format");
			
			x0 = Double.parseDouble(array[0]);
			xf = Double.parseDouble(array[1]);
			y0 = Double.parseDouble(array[2]);
			yf = Double.parseDouble(array[3]);
			phasePlot.repaint();
		} catch (IllegalArgumentException e) {
			JOptionPane.showMessageDialog(phasePlot, e.getMessage(), e.getClass().getSimpleName(),
					JOptionPane.ERROR_MESSAGE);
		} catch (Exception e) {}
	}
	
	// the method that prints the phase state of the system
	private void printPhaseVariables ()
	{
		// execute in a worker thread
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {
			@Override
			public void run () {
				try {
					// launch the file dialog for the file name
					FileDialog dialog = new FileDialog (PhasePlotter.this);
					dialog.setVisible(true);
					
					// get the chosen file
					String fileName = dialog.getFile();
					
					// check that a file has been chosen
					if (fileName != null && !fileName.equals("")) {
						// create the file and instantiate a print writer
						java.io.File file = new java.io.File(
								dialog.getDirectory(), fileName);
						if (!file.exists())
							file.createNewFile();
						java.io.FileWriter writer = new java.io.FileWriter(file);
						
						// prepare the parameters
						DecimalFormat df = new DecimalFormat ("#.######");
						String paramStr = "";
						for (int i = 0; i < par.length; ++i)
							paramStr += par [i] + " = " + 
									df.format(var[parIndices[i]]) + ", ";
						paramStr = paramStr.substring(0, paramStr.lastIndexOf(", "));
						paramStr += "\r\n";
						
						// show the variable order
						for (String xx : x)
							paramStr += xx + "\t";
						
						// write the paramStr to file
						writer.write(paramStr + "\r\n");
						
						// get the phase variables and write to the chosen file
						for (int i = 0; i < xVal.length; ++i){
							for (int j = 0; j < xVal[i].length; ++j)
								writer.write (df.format(xVal[i][j]) + "\t");
							
							writer.write("\n");
						}
						
						// show success message
						javax.swing.JOptionPane.showMessageDialog(PhasePlotter.this, 
								"The phase plot has been printed successfully",
								"PRINT DIALOG", javax.swing.JOptionPane.ERROR_MESSAGE);
						
						writer.flush();
						writer.close();
					}
				} catch (java.io.IOException e) {
					javax.swing.JOptionPane.showMessageDialog(PhasePlotter.this, 
							e.getMessage(), e.getClass().getSimpleName(),
							javax.swing.JOptionPane.ERROR_MESSAGE);
				}
			}
		});
		executor.shutdown();
	}	// end method print Phase variables
	
	// method to show the parameter on the paramArea
	private void showParameters ()
	{
		paramArea.setText("");
		DecimalFormat df = new DecimalFormat ("#.####");
		
		for (int i = 0; i < par.length; ++i)
			paramArea.append(par[i] + " = " + df.format(
					val[parIndices[i]]) + "   ");
		paramArea.append("\n");
		for (int i = 0; i < x.length; ++i)
			paramArea.append(x[i] + " = " + df.format(
					xVal[plotCounter - 1][xIndices[i]]) + "   ");
		
		// show the analysis variables
		paramArea.append("\nTime: " + df.format(
				(dt * (plotCounter - steadyStateCounter))));
	}
    
    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
    	try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        } catch (InstantiationException e) {
        	e.printStackTrace();
        } catch (IllegalAccessException e) {
        	e.printStackTrace();
        } catch (javax.swing.UnsupportedLookAndFeelException e) {
        	e.printStackTrace();
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
            	new PhasePlotter ().setVisible(true);
            }
        });
    }

    // Frame constructor
    public PhasePlotter ()
    {
    	super ("Ph.D models");
    	
    	// initialize the plot variables
    	initializePlotter ();
    	
    	// create the popup menu
    	popupMenu = new javax.swing.JPopupMenu("Ph.D model popup menu");
    	
    	steadyMenu = new javax.swing.JCheckBoxMenuItem("Plot Steady State", showSteadyState);
    	iteratedMenu = new javax.swing.JCheckBoxMenuItem("Use Iterated Plot", iteratedPlot);
    	javax.swing.JMenuItem createMenu = new javax.swing.JMenuItem("Create New Model");
    	javax.swing.JMenuItem selectMenu = new javax.swing.JMenuItem("Select Model from List");
    	javax.swing.JMenuItem nMenu = new javax.swing.JMenuItem("Change n iterations");
    	javax.swing.JMenuItem plotVarMenu = new javax.swing.JMenuItem("Change Plot variables");
    	javax.swing.JMenuItem varRangeMenu = new javax.swing.JMenuItem("Change variable range");
    	javax.swing.JMenuItem printMenu = new javax.swing.JMenuItem("Print Phase variables");
    	javax.swing.JMenuItem setParamMenu = new javax.swing.JMenuItem("Set Parameter values");
    	javax.swing.JMenuItem paramIncrementMenu = new javax.swing.JMenuItem(
    			"Set Parameter Increment");
    	javax.swing.JMenuItem changeParamMenu = new javax.swing.JMenuItem("Change Control Parameter");
    	
    	
    	popupMenu.add(new javax.swing.JSeparator());
    	popupMenu.add(createMenu);
    	popupMenu.add(selectMenu);
    	popupMenu.add(new javax.swing.JSeparator());
    	popupMenu.add(plotVarMenu);
    	popupMenu.add(varRangeMenu);
    	popupMenu.add(new javax.swing.JSeparator());
    	popupMenu.add(setParamMenu);
    	popupMenu.add(changeParamMenu);
    	popupMenu.add(paramIncrementMenu);
    	popupMenu.add(new javax.swing.JSeparator());
    	popupMenu.add(steadyMenu);
    	popupMenu.add(iteratedMenu);
    	popupMenu.add(nMenu);
    	popupMenu.add(new javax.swing.JSeparator());
    	popupMenu.add(printMenu);
    	popupMenu.add(new javax.swing.JSeparator());
    	
    	createMenu.addActionListener(new java.awt.event.ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				// TODO Auto-generated method stub
				showModelDialog ();
			}
    		
    	});
    	
    	selectMenu.addActionListener(new java.awt.event.ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				// TODO Auto-generated method stub
				selectModel ();
			}
    		
    	});
    	
    	setParamMenu.addActionListener(new java.awt.event.ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				// TODO Auto-generated method stub
				setParameter ();
			}
    		
    	});
    	
    	changeParamMenu.addActionListener(new java.awt.event.ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				// TODO Auto-generated method stub
				changeParameter ();
			}
    		
    	});
    	
    	paramIncrementMenu.addActionListener(new java.awt.event.ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				// TODO Auto-generated method stub
				setParameterIncrement ();
			}
    		
    	});
    	
    	steadyMenu.addActionListener (new java.awt.event.ActionListener() {
			@Override
			public void actionPerformed (ActionEvent e) {
				// TODO Auto-generated method stub
				setSteadyState ();
			}
    	});

    	iteratedMenu.addActionListener (new java.awt.event.ActionListener() {
			@Override
			public void actionPerformed (ActionEvent e) {
				// TODO Auto-generated method stub
				setIteratedPlot ();
			}
    	});
    	
    	nMenu.addActionListener(new java.awt.event.ActionListener () {
    		public void actionPerformed (java.awt.event.ActionEvent e) {
    			// prompt for input
    			try {
    				n = Integer.parseInt(JOptionPane.showInputDialog(PhasePlotter.this, 
    						"Enter the new number of iterations"));
    				xVal = new double [n] [xIndices.length];
    				sx = new double [n];
    				sy = new double [n];
    				getSolution ();
    			} catch (NumberFormatException nfe) {
					JOptionPane.showMessageDialog(PhasePlotter.this, nfe.getMessage(),
							nfe.getClass().getSimpleName(), JOptionPane.ERROR_MESSAGE);
				}
    		}
    	});

    	plotVarMenu.addActionListener(new java.awt.event.ActionListener () {
    		public void actionPerformed (java.awt.event.ActionEvent e) {
    			changePlotVariables ();
    		}
    	});
    	
    	varRangeMenu.addActionListener(new java.awt.event.ActionListener () {
    		public void actionPerformed (java.awt.event.ActionEvent e) {
    			changeVariableRange();
    		}
    	});
    	
    	printMenu.addActionListener(new java.awt.event.ActionListener () {
    		public void actionPerformed (java.awt.event.ActionEvent e) {
    			printPhaseVariables();
    		}
    	});
    	
    	phasePlot = new javax.swing.JPanel() {
			private static final long serialVersionUID = 1L;

			// define this paintComponent to a method
    		public void paintComponent (Graphics g)
    	 	{
    	 		if (!multiplePlot)
    	 			super.paintComponent(g);
    	 		
		 		// plot the phase state -- checking for special plot
    	 		g.setColor(Color.black);
    	 		if (specialPlot) {
    	 			for (int i = steadyStateCounter; i < plotCounter; ++i)
        	 			g.drawLine(toPx(sx[i - 1]), toPy(sy[i - 1]), 
        	 					toPx(sx[i]), toPy(sy[i]));
    	 		} else {
    	 			for (int i = steadyStateCounter; i < plotCounter; ++i) {
        	 			g.drawLine(
        	 					toPx(xVal[i - 1][plotIndices[0]]), 
        	 					toPy(xVal[i - 1][plotIndices[1]]),
        	 					toPx(xVal[i][plotIndices[0]]), 
        	 					toPy(xVal[i][plotIndices[1]]));
        	 		}
    	 		}
    	 		
    	 		// draw the inner rectangle
		 		g.drawRect(xStart, yStart, width, height);
    		 	
		 		// show the legend
				g.drawString("" + yf, xStart - 25, yStart + 5);
				g.drawString("" + y0, xStart - 25, yStart + height + 5);
				g.drawString("" + x0, xStart - 5, yStart + height + 15);
				g.drawString("" + xf, xStart + width - 10, yStart + height + 15);
				
				// show the plot variables
				String xAxis = (specialPlot) ? specialLegendX : plotVar[0];
				String yAxis = (specialPlot) ? specialLegendY : plotVar[1];
				g.drawString(xAxis, (int)(xStart + 0.5 * width), yStart + height + 20);
				g.drawString(yAxis, 15, (int)(yStart + 0.5 * height));
				
				// show minor marks
				for (int i = 0; i < 21; ++i) {
					// divide y-axis
					g.drawString("-", xStart, yStart + (int) (i * (height) / 20) + 4);
					g.drawString("-", xStart + width - 2, yStart + (int) (i * (height) / 20) + 4);
					g.drawString("'", xStart + (int) (i * width / 20) - 1, yStart + 10);
					g.drawString("'", xStart + (int) (i * width / 20) - 1, yStart + (height) + 7);
				}
    	 	}
    	};
    	
    	jPanel2 = new javax.swing.JPanel();
		paramArea = new javax.swing.JTextArea();
		paramArea.setEditable(false);
		paramArea.setWrapStyleWord(true);
		javax.swing.JScrollPane paramPane = new javax.swing.JScrollPane(paramArea);
		
        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        phasePlot.setBackground(new java.awt.Color(255, 255, 255));
        phasePlot.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        phasePlot.setPreferredSize(new java.awt.Dimension(500, 500));

        phasePlot.addMouseListener(new java.awt.event.MouseAdapter () {

			public void mouseClicked(MouseEvent e) {
				// TODO Auto-generated method stub
				// get the clicked point
				double xx = toX(e.getX());
				double yy = toY(e.getY());
				if (xx >= x0 && xx <= xf && yy >= y0 && yy <= yf)
					setPlotVariableInitial (toX(e.getX()), toY(e.getY()));
			}

			public void mouseReleased(MouseEvent e) {
				// TODO Auto-generated method stub
				if (e.isPopupTrigger())
					popupMenu.show(e.getComponent(), e.getX(), e.getY());
			}
        });
        
        // respond to keyboard actions
        paramArea.addKeyListener(new java.awt.event.KeyAdapter () {
			public void keyTyped(KeyEvent e) {
				// TODO Auto-generated method stub
				switch (e.getKeyChar()) {
					case 'i': case 'I':
						// prompt for initial condition
						setVariablesInitial();
						break; 
					case 'x' : case 'X':
						setParameter();
						break;
					case 'c' : case 'C':
						if (!backward && iteratedPlot && plotCounter < n)
							++plotCounter;
						else if (backward && iteratedPlot && plotCounter > 0)
							--plotCounter;
						repaint();
						showParameters();
						break;
					case 'v' : case 'V':
						if (!backward && iteratedPlot && plotCounter < n - 10)
							plotCounter += 10;
						else if (backward && iteratedPlot && plotCounter > 10)
							plotCounter -= 10;
						repaint();
						showParameters();
						break;
					case 'b' : case 'B':
						if (!backward && iteratedPlot && plotCounter < n - 100)
							plotCounter += 100;
						else if (backward && iteratedPlot && plotCounter > 100)
							plotCounter -= 100;
						repaint();
						showParameters();
						break;
					case 'n' : case 'N':
						if (!backward && iteratedPlot && plotCounter < n - 1000)
							plotCounter += 1000;
						else if (backward && iteratedPlot && plotCounter > 1000)
							plotCounter -= 1000;
						repaint();
						showParameters();
						break;
					case 'm' : case 'M':
						if (!backward && iteratedPlot && plotCounter < n - 10000)
							plotCounter += 10000;
						else if (backward && iteratedPlot && plotCounter > 10000)
							plotCounter -= 10000;
						repaint();
						showParameters();
						break;
					case '-':
						backward = true;
						if (plotCounter > steadyStateCounter)
							--plotCounter;
						repaint();
						showParameters();
						break;
					case '+': case '=':
						backward = false;
						if (plotCounter < n)
							++plotCounter;
						repaint();
						break;
					case '1' : case '!':
						if (!backward)
							val [controlParam[0]] += dp;
						else
							val [controlParam[0]] -= dp;
						
						plotCounter = steadyStateCounter;
						getSolution();
						showParameters();
						break;
					case '2' : case '@':
						if (!backward)
							val [controlParam[1]] += dp;
						else
							val [controlParam[1]] -= dp;
						
						plotCounter = steadyStateCounter;
						getSolution();
						showParameters();
						break;
					case 'p' : case 'P':
						setParameterIncrement ();
						showParameters();
						break;
					default:
				}
			}
		});
        
        javax.swing.GroupLayout jPanel2Layout = new javax.swing.GroupLayout(jPanel2);
        jPanel2.setLayout(jPanel2Layout);
        jPanel2Layout.setHorizontalGroup(jPanel2Layout.createSequentialGroup()
            .addComponent(paramPane, 0, 500, Short.MAX_VALUE)
        );
        jPanel2Layout.setVerticalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addComponent(paramPane, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE)
                )
        );
        
        javax.swing.GroupLayout phasePlotLayout = new javax.swing.GroupLayout(phasePlot);
        phasePlot.setLayout(phasePlotLayout);
        phasePlotLayout.setHorizontalGroup(
            phasePlotLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 500, Short.MAX_VALUE)
        );
        phasePlotLayout.setVerticalGroup(
            phasePlotLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 500, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                		.addContainerGap()
                        .addComponent(jPanel2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(phasePlot, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(5, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(phasePlot, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(10, 10, 10)
                .addComponent(jPanel2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(0, 0, 0))
        );
        
        pack();
        //setLocationRelativeTo (null);
        setResizable (false);
    	
	}	// end Frame's constructor
    
    /**
	 * This method converts x-value to the corresponding pixels
	 * @param x
	 * @return
	 */
	private int toPx (double x)
	{
		return (int) (xStart + (x - x0) / (xf - x0) * width);
	}
	
	/** 
	 * This method converts a y-value to equivalent y on the screen
	 * @param y
	 * @return
	 */
	private int toPy (double y)
	{
		return (int) (yStart + (y - yf) / (y0 - yf) * height);
	}
	
	private double toX (double Px)
	{
		return (x0 + (xf - x0) * (Px - xStart) / width);
	}
	
	private double toY (double Py)
	{
		return (yf + (y0 - yf) * (Py - yStart) / height);
	}
	
	// inner class
	private class PlotModel {
		public String modelName;
		public ExpNode [] fx;
		public String [] x;
		public String [] var;
		public double [] val;
		public int [] coord;
		
		public PlotModel (String m, ExpNode [] f, String [] x, String [] var,
				double [] val, int [] c)
		{
			this.modelName = m;
			this.fx = f;
			this.x = x;
			this.var = var;
			this.val = val;
			this.coord = c;
		}
		
		@Override
		public boolean equals (Object obj)
		{
			// cast obj to model
			if (obj instanceof PlotModel) {
				PlotModel model = (PlotModel) obj;
				if (model.modelName.equals(modelName))
					return true;
			}
			return false;
		}
		
		@Override
		public String toString () {
			return modelName;
		}
	}
	
	
    // Variables declaration - do not modify 
    private javax.swing.JPanel phasePlot;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPopupMenu popupMenu;
    private javax.swing.JTextArea paramArea;
    private javax.swing.JCheckBoxMenuItem steadyMenu;
    private javax.swing.JCheckBoxMenuItem iteratedMenu;
    
    // End of variables declaration 

}
