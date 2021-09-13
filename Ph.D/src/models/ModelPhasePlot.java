package models;

import java.awt.Color;
import java.awt.FileDialog;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.MouseEvent;
import java.text.DecimalFormat;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.JOptionPane;

public class ModelPhasePlot extends javax.swing.JFrame
{
	private static final long serialVersionUID = 1L;
	
	// instance variables
	private static int n = 100000;
	private double dt = 0.01;
	
	private double p10 = 0.0, p20 = 0.0, p1f = 2.0, p2f = 2.0;
	
	private double x0 = -3.0;
	private double xf = 3.0;
	private double y0 = -3.0;
	private double yf = 3.0;
	
	private static int xStart = 50;
	private static int yStart = 50;
	private static int width = 400;
	private static int height = 400;
	
	// phase plot variables and parameter definition
	private String [] variableList;
	private String [] parameterList;
	private int xIndex = 0;
	private int yIndex = 1;
	private int p1Index = 0;
	private int p2Index = 1;
	private double [] param;
	private double [][] v;
	
	private int model = -1;
	private String modelList = "";
	
	private boolean multiplePlot = false;
	private boolean showSteadyState = false;
	private int steadyStateIterate = 1;
	
	private boolean specialPlot = false;
	private String specialX = "", specialY = "";
	private String specialLegendX = "", specialLegendY = "";
	private double [] sx, sy;
	
    // initialization method
	private void initializePhasePlot ()
	{
		// continuously request the user for model
		while (true) {
			try {
				modelList = "Choose \n" +
						"1 for Model 17\n" +
						"2 for Model 31\n" +
						"3 for Sekikawa2010\n" +
						"4 for Ueta2003\n" +
						"5 for Ueta2004\n" +
						"6 for Original van der Pol oscillator or \n" +
						"-1 to quit\n";
				String input = javax.swing.JOptionPane.showInputDialog (null, modelList,
						"Model Type", javax.swing.JOptionPane.INFORMATION_MESSAGE);
				model = Integer.parseInt(input);
				choseModel (model);
				break;
			} catch (NumberFormatException e) {
				javax.swing.JOptionPane.showMessageDialog(null,
						"Illegal input! Press -1 to quit",
						e.getClass().getSimpleName(), 
						javax.swing.JOptionPane.ERROR_MESSAGE);
			} catch (IllegalArgumentException e) {
				javax.swing.JOptionPane.showMessageDialog(null, e.getMessage(), 
						e.getClass().getSimpleName(), 
						javax.swing.JOptionPane.ERROR_MESSAGE);
			} 
		}
		
	}
	
	// method that allows the user to choose the model of choice
	private void choseModel (int model) throws IllegalArgumentException
	{
		switch (model) 
		{
			case 1:
				// model 17
				setTitle ("Phase portrait for Model 17");
				// variable declaration
				variableList = new String [] {"x1", "y1", "x2", "y2", "phi", "t"};
				parameterList = new String [] {"d", "g1", "g2", "k1", "k2", "s", "e3"};
				
				// instantiate variables v and parameters
				v = new double [n] [variableList.length];
				v[0] = new double [] {0.6, 0.4, 0.7, 0.3, 0.01, 0.0};
				break;
				
			case 2:
				// model 31
				setTitle ("Phase portrait for Model 31");
				// variable declaration
				variableList = new String [] {"x1", "y1", "x2", "y2", "phi", "t"};
				parameterList = new String [] {"g1", "g2", "k1", "k2", "R1", "R2", "e3", "n"};
				
				// instantiate variables v and parameters
				v = new double [n] [variableList.length];
				v[0] = new double [] {0.6, 0.4, 0.7, 0.3, 0.01, 0.0};
				break;
				
			case 3:
				// seikikawa2010
				setTitle ("Phase portrait for Seikikawa 2010");
				// variable declaration
				variableList = new String [] {"x", "y", "z", "t"};
				parameterList = new String [] {"e", "k1", "k2", "k3", "B0", "B1"};
				
				// instantiate variables v and parameters
				v = new double [n] [variableList.length];
				v[0] = new double [] {0.6, 0.4, 0.7, 0.0};
				break;
			case 4:
				setTitle ("Phase portrait for Ueta 2003");
				// assymetrically coupled bvp oscillators {ueta2003}
				variableList = new String [] {"x1", "y1", "x2", "y2", "t"};
				parameterList = new String [] {"g", "k1", "k2", "d"};
				
				// instantiate variables v and parameters
				v = new double [n] [variableList.length];
				v[0] = new double [] {0.5, 0.5, 0.5, 0.5, 0.0};
				break;
				
			case 5:
				// Definition for a resistively current coupled BVP	ueta2004
				setTitle ("Phase portrait for Ueta2004");
				// define variableList and parameterList
				variableList = new String[] {"x1", "y1", "x2", "y2", "t"};
				parameterList = new String[] {"g", "k", "d"};
				
				// instantiate variables v and parameters
				v = new double [n] [variableList.length];
				v[0] = new double [] {0.5, 0.5, 0.5, 0.5, 0.0};
				break;
				
			case 6:
				setTitle ("Phase portrait for Original van der Pol Oscillator");
				// definition for the original van der Pol equation
				variableList = new String [] {"x", "y", "t"};
				parameterList = new String [] {"c", " "};
				
				// instantiate variables v and parameters
				v = new double [n] [variableList.length];
				v[0] = new double [] {0.5, 0.5, 0.0};
				break;
			
			case -1:
				System.exit (0);
				break;
				
			default:
				throw new IllegalArgumentException ("The chosen model does not exit!");
		}
		
		// set the other variables
		param = new double [parameterList.length];
		
		// assign special variables
		sx = new double [v.length];
		sy = new double [v.length];
	}
	
	// method to set the variable initial condition
	private void setVariableInitialCondition (double x, double y)
	{
		v[0][xIndex] = x;
		v[0][yIndex] = y;
		getSolution ();
	}	// end method setVariableInitialCondition

    
    // method to solve for the system and plot 
    private void getSolution ()
    {
    	// iterate n-times
    	for (int i = 1; i < v.length; ++i)
    		v[i] = rungeKutta (v[i - 1]);
    	
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
        		int index1 = Integer.parseInt(specialX.substring(0, specialX.indexOf("+")));
        		int index2 = Integer.parseInt(specialX.substring(
        				specialX.indexOf("+"), specialX.length()));
        		int indexY = Integer.parseInt(specialY);
        		
        		// check the entered indices
        		int max = variableList.length;
        		if (index1 >= max || index2 >= max || indexY >= max)
        			throw new IllegalArgumentException ("The entered index is out of range");
        		
        		// set the values for the special plot
        		for (int i = 0; i < sx.length; ++i) {
        			sx[i] = v[i][index1] + v[i][index2];
        			sy[i] = v[i][indexY];
        		}
        		// set the legend for the special axis
        		specialLegendX = variableList[index1] + " + " + variableList[index2];
        		specialLegendY = variableList[indexY];
        	}
        	
        	// assign for y
        	if (specialY.contains("+")) {
        		int index1 = Integer.parseInt(specialY.substring(0, specialY.indexOf("+")));
        		int index2 = Integer.parseInt(specialY.substring(
        				specialY.indexOf("+"), specialY.length()));
        		int indexX = Integer.parseInt(specialX);
        		
        		// check the entered indices
        		int max = variableList.length;
        		if (index1 >= max || index2 >= max || indexX >= max)
        			throw new IllegalArgumentException ("The entered index is out of range");
        		
        		
        		// set the values for the special plot
        		for (int i = 0; i < sy.length; ++i) {
        			sy[i] = v[i][index1] + v[i][index2];
        			sx[i] = v[i][indexX];
        		}
        		
        		// set the legend for the special axis
        		specialLegendY = variableList[index1] + " + " + variableList[index2];
        		specialLegendX = variableList[indexX];
        	}
        	
    	} catch (IllegalArgumentException e) {
			javax.swing.JOptionPane.showMessageDialog(phasePlot, e.getMessage(),
					"INPUT ERROR", javax.swing.JOptionPane.ERROR_MESSAGE);
		}
    }
    
    // Runge Kutta fourth order numerical integeration
  	private double [] rungeKutta (double [] v)
  	{
  		int l = v.length;
  		double [] c1, c2, c3, c4;
  		
  		// initialize the intermediate steps
  		c1 = new double [l];
  		c2 = new double [l];
  		c3 = new double [l];
  		c4 = new double [l];
  		
  		c1 = f(v);
  		
  		for (int i = 0; i < l; ++i)
  			c2[i] = v[i] + dt * c1[i] / 2;
  		c2 = f(c2);
  		
  		for (int i = 0; i < l; ++i)
  			c3[i] = v[i] + dt * c2[i] / 2;
  		c3 = f(c3);
  		
  		for (int i = 0; i < l; ++i)
  			c4[i] = v[i] + dt * c3[i];
  		c4 = f(c4);
  		
  		for (int i = 0; i < l; ++i)
  			c1[i] = v[i] + dt * (c1[i] + 2 * (c2[i] + c3[i]) + c4[i]) / 6.0;
  		
  		return c1;
  	}
 	
 	/**
 	 * This method defines the equation, set the variables and define the parameters
 	 * for the system.
 	 * @param g
 	 * @param k
 	 * @param v
 	 * @return
 	 */
 	public double [] f (double [] v)
 	{
 		int l = v.length;
 		double [] vv = new double [l];
 		double e1 = -0.00005;
 		
		switch (model) 
 		{
 		case 1:	// model 17
 			vv[0] = -v[1] + Math.tanh(param[1] * v[0]) - 
 					param[0] * (e1 + param[6] * v[4] * v[4]) * (v[0] - v[2]);
 			vv[1] = v[0] - param[3] * v[1];
 			vv[2] = -v[3] + Math.tanh(param[2] * v[2]) - 
 					param[0] * (e1 + param[6] * v[4] * v[4]) * (v[2] - v[0]);
 			vv[3] = v[2] - param[4] * v[3];
 			vv[4] = param[5] * (v[0] - v[2]);
 			vv[5] = 1.0;
 			break;
 			
 		case 2: // model 31
 			double ff = 1.0 * (param[4] * v[1] - param[5] * v[3]) /
 	 				(1.0 + (e1 + param[6] * v[4] * v[4]) * (param[4] + param[5]));
 	 		
 	 		vv[0] = -v[1] + Math.tanh(param[0] * v[0]);
 	 		vv[1] = v[0] - param[2] * v[1] + param[2] * ff * (e1 + param[6] * v[4] * v[4]);
 	 		vv[2] = -v[3] + Math.tanh(param[1] * v[2]);
 	 		vv[3] = v[2] - param[3] * v[3] + param[3] * -ff * (e1 + param[6] * v[4] * v[4]);
 	 		vv[4] = param[7] * ff;
 	 		vv[5] = 1.0;
 	 		break;
 			
 		case 3:
 			vv[0] = (v[0] * (1 - v[0] * v[0]) + v[1] + v[2]) / param[0];
 			vv[1] = -v[0] - param[1] * v[1] + param[4];
 			vv[2] = param[3] * (-v[0] - param[2] * v[2] + param[5]);
 			vv[3] = 1.0;
 			break;
 			
 		case 4:
 			double n = 1.0 / (1 + param[3] * param[1]);
 	 		
 	 		vv[0] = -v[1] + Math.tanh(param[0] * v[0]);
 	 		vv[1] = v[0] - param[1] * v[1] + param[3] * param[1] * n * (
 	 				param[1] * v[1] - v[2]);
 	 		vv[2] = -v[3] + Math.tanh(param[0] * v[2]) + param[3] * n * (
 	 				param[1] * v[1] - v[2]);
 	 		vv[3] = v[2] - param[2] * v[3];
 	 		vv[4] = 1.0;
 	 		break;
 			
 		case 5:
 			vv[0] = -v[1] + Math.tanh(param[0] * v[0]);
 	 		vv[1] = v[0] - param[1] * v[1] + param[2] * param[1] * (v[1] - v[3]);
 	 		vv[2] = -v[3] + Math.tanh(param[0] * v[2]);
 	 		vv[3] = v[2] - param[1] * v[3] + param[2] * param[1] * (v[3] - v[1]);
 	 		vv[4] = 1.0;
 			break;
 			
 		case 6:
 			vv[0] = v[1];
 			vv[1] = -v[0] - param[0] * (v[0] * v[0] - 1) * v[1];
 			vv[2] = 1;
 			break;
 			
		default:
			throw new IllegalArgumentException ("Illegal chosen model.");
 		}
		
 		return vv;
 	}
 	
	// method to change the plot variables
	private void changePlotVariables ()
	{
		// set the special plot to false
		specialPlot = false;
		
		try {
			// prepare the output message
			String varStr = "";
			for (int i = 0; i < variableList.length; ++i)
				varStr += i + "   is for   " + variableList[i] + "\n";
			
			// obtain user's input
			String input = javax.swing.JOptionPane.showInputDialog(phasePlot, 
					"Enter the new variables as 0,1 where \n" + varStr
					+ " and variable 0 is against variable 1",
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
				int index1 = Integer.parseInt(str[0]);
				int index2 = Integer.parseInt(str[1]);
				
				// confirm that indices are valid
				if (index1 >= variableList.length)
					throw new IllegalArgumentException ("The first chosen index is invalid");
				if (index2 >= variableList.length)
					throw new IllegalArgumentException ("The second chosen index is invalid");
				
				// assign values
				xIndex = index2;
				yIndex = index1;
			}
			
			phasePlot.repaint();
		} catch (IllegalArgumentException e) {
			javax.swing.JOptionPane.showMessageDialog(phasePlot, e.getMessage(),
					"INPUT ERROR", javax.swing.JOptionPane.ERROR_MESSAGE);
		} catch (Exception e) {}
	}
	
	// method to change the control parameters
	private void changeControlParameters ()
	{
		try {
			// prepare the output message
			String paramStr = "";
			for (int i = 0; i < parameterList.length; ++i)
				paramStr += i + "   is for   " + parameterList[i] + "\n";
			
			// obtain user's input
			String input = javax.swing.JOptionPane.showInputDialog(phasePlot, 
					"Enter the new parameters as 0,1 where \n" + paramStr
					+ " and parameter 0 is against variable 1",
					javax.swing.JOptionPane.YES_NO_OPTION);
			
			// split input around comma
			String [] str = input.split(",");
			
			if (str.length != 2)
				throw new IllegalArgumentException ("Invalid input format!!!");
			
			// convert input
			int index1 = Integer.parseInt(str[0]);
			int index2 = Integer.parseInt(str[1]);

			// confirm that indices are valid
			if (index1 >= parameterList.length && index2 >= parameterList.length)
				throw new IllegalArgumentException ("Chosen index is invalid");
			
			// assign values
			p1Index = index1;
			p2Index = index2;
			
			setParameterChange();
		} catch (IllegalArgumentException e) {
			javax.swing.JOptionPane.showMessageDialog(phasePlot, e.getMessage(),
					"INPUT ERROR", javax.swing.JOptionPane.ERROR_MESSAGE);
		} catch (Exception e) {}
	}
	
	// method to set the parameter values
	private void setParameterValues ()
	{
		try {
			String format = "";
			for (String str : parameterList)
				format += str + ",";
			format = format.substring(0, format.length() - 1);
			
			// prompt the user for values
			String input = JOptionPane.showInputDialog(phasePlot, 
					"Enter the parameter values as: \n" + format, JOptionPane.YES_NO_OPTION);
			
			// split input around commas
			String[] values = input.split(",");
			if (values.length != parameterList.length)
				throw new IllegalArgumentException ("Invalid Entry! The number of entered values " +
						"is not equal to parameter list");
			for (int i = 0; i < values.length; ++i)
				param[i] = Double.parseDouble(values[i]);
			
			// update the values on the sliders
			setParameterChange();
		} catch (IllegalArgumentException e) {
			javax.swing.JOptionPane.showMessageDialog(phasePlot, e.getMessage(),
					"INPUT ERROR", javax.swing.JOptionPane.ERROR_MESSAGE);
		} catch (Exception e) {}
	}
	
	// method to change the display on the slider during parameter change
	private void setParameterChange ()
	{
		// show the current parameter values and label
		p1Label.setText(parameterList[p1Index]);
		p2Label.setText(parameterList[p2Index]);
		p1ValueLabel.setText(new DecimalFormat ("#.####").format (param[p1Index]));
		p2ValueLabel.setText(new DecimalFormat ("#.####").format (param[p2Index]));
		
		// change value of the slider
		p1Slider.setValue((int) (p1Slider.getMinimum() + 
				(param[p1Index] - p10) * p1Slider.getMaximum() / (p1f - p10)));
		p2Slider.setValue((int) (p2Slider.getMinimum() + 
				(param[p2Index] - p10) * p1Slider.getMaximum() / (p2f - p20)));
	}
	
	// method to change the parameter range
	private void changeParameterRange ()
	{
		try {
			String input = JOptionPane.showInputDialog(phasePlot, "Enter paramter range for p1 and p2 " +
					"as P1o,P1f,p2o,p2f format", "Parameter Range Alteration Dialog", 
					JOptionPane.INFORMATION_MESSAGE);
			
			// split string around commas into parts
			String [] array = input.split(",");
			
			// confirm proper entry
			if (array.length != 4)
				throw new IllegalArgumentException ("The entered input does not " +
						"confirm to P1o,P1f,p2o,p2f format");
			
			p10 = Double.parseDouble(array[0]);
			p1f = Double.parseDouble(array[1]);
			p20 = Double.parseDouble(array[2]);
			p2f = Double.parseDouble(array[3]);
			evaluateP1();
			evaluateP2();
		} catch (IllegalArgumentException e) {
			JOptionPane.showMessageDialog(phasePlot, e.getMessage(), e.getClass().getSimpleName(),
				JOptionPane.ERROR_MESSAGE);
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
					FileDialog dialog = new FileDialog (ModelPhasePlot.this);
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
						for (int i = 0; i < parameterList.length; ++i)
							paramStr += parameterList [i] + " = " + 
									df.format(param[i]) + ", ";
						paramStr = paramStr.substring(0, paramStr.lastIndexOf(", "));
						paramStr += "\r\n";
						
						// show the variable order
						for (String var : variableList)
							paramStr += var + "\t";
						
						// write the paramStr to file
						writer.write(paramStr + "\r\n");
						
						// get the phase variables and write to the chosen file
						for (int i = 0; i < v.length; ++i){
							for (int j = 0; j < v[i].length; ++j)
								writer.write (df.format(v[i][j]) + "\t");
							
							writer.write("\n");
						}
						
						// show success message
						javax.swing.JOptionPane.showMessageDialog(ModelPhasePlot.this, 
								"The phase plot has been printed successfully",
								"PRINT DIALOG", javax.swing.JOptionPane.ERROR_MESSAGE);
						
						writer.flush();
						writer.close();
					}
				} catch (java.io.IOException e) {
					javax.swing.JOptionPane.showMessageDialog(ModelPhasePlot.this, 
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
		for (int i = 0; i < parameterList.length; ++i)
			paramArea.append(parameterList[i] + " = " + df.format(param[i]) + "   ");
		paramArea.append("\n");
		for (int i = 0; i < variableList.length; ++i)
			paramArea.append(variableList[i] + "=" + df.format(v[0][i]) + "   ");
	}
	
	// method to change param minor tics
	private void changeParamMinorTick()
	{
		
	}
	
	// method to change parameter major ticks
	private void changeParamMajorTick()
	{
		
	}
	
	// method to change pa
	private void evaluateP2 ()
	{
        // TODO add your handling code here:
    	// obtain the values of k- and g-
    	param[p2Index] = p20 + (p2f - p20) * p2Slider.getValue() / p2Slider.getMaximum();
    	p2ValueLabel.setText(new DecimalFormat ("#.####").format (param[p2Index]));
    	getSolution();
    	//initialCondition ();
    }                                    

    private void evaluateP1 ()
    {
        // TODO add your handling code here:
    	param[p1Index] = p10 + (p1f - p10) * p1Slider.getValue() / p1Slider.getMaximum();
    	p1ValueLabel.setText(new DecimalFormat ("#.####").format(param[p1Index]));
    	getSolution();
    	//initialCondition ();
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
            	new ModelPhasePlot ().setVisible(true);
            }
        });
    }

    // Frame constructor
    public ModelPhasePlot ()
    {
    	super ("Ph.D models");
    	
    	// create the popup menu
    	popupMenu = new javax.swing.JPopupMenu("Ph.D model popup menu");
    	
    	final javax.swing.JCheckBoxMenuItem steadyMenu = 
    			new javax.swing.JCheckBoxMenuItem("Plot Steady State", showSteadyState);
    	javax.swing.JMenuItem nMenu = new javax.swing.JMenuItem("Change n iterations");
    	javax.swing.JMenuItem plotVarMenu = new javax.swing.JMenuItem("Change Plot variables");
    	javax.swing.JMenuItem paramMenu = new javax.swing.JMenuItem("Change control parameters");
    	javax.swing.JMenuItem varRangeMenu = new javax.swing.JMenuItem("Change variable range");
    	javax.swing.JMenuItem paramRangeMenu = new javax.swing.JMenuItem("Change parameter range");
    	javax.swing.JMenuItem setParamMenu = new javax.swing.JMenuItem("Set Parameter values");
    	javax.swing.JMenuItem paramMinorMenu = new javax.swing.JMenuItem("Change parameter minor ticks");
    	javax.swing.JMenuItem paramMajorMenu = new javax.swing.JMenuItem("Change parameter major ticks");
    	javax.swing.JMenuItem printMenu = new javax.swing.JMenuItem("Print Phase variables");
    	
    	popupMenu.add(new javax.swing.JSeparator());
    	popupMenu.add(plotVarMenu);
    	popupMenu.add(varRangeMenu);
    	popupMenu.add(new javax.swing.JSeparator());
    	popupMenu.add(paramMenu);
    	popupMenu.add(paramRangeMenu);
    	popupMenu.add(setParamMenu);
    	popupMenu.add(new javax.swing.JSeparator());
    	popupMenu.add(paramMinorMenu);
    	popupMenu.add(paramMajorMenu);
    	popupMenu.add(new javax.swing.JSeparator());
    	popupMenu.add(steadyMenu);
    	popupMenu.add(nMenu);
    	popupMenu.add(printMenu);
    	popupMenu.add(new javax.swing.JSeparator());
    	
    	steadyMenu.addActionListener (new java.awt.event.ActionListener() {
			@Override
			public void actionPerformed (ActionEvent e) {
				// TODO Auto-generated method stub
				
				// toggle showSteadyState
				showSteadyState = !showSteadyState;
				
				// prompt the user for entry
				if (showSteadyState) {
					try {
						steadyStateIterate = Integer.parseInt(JOptionPane.showInputDialog(
								ModelPhasePlot.this, "Enter the steady state iterate: "));
						
						// check for valid input for the steady state
						if (steadyStateIterate > n || steadyStateIterate < 1)
							throw new IllegalArgumentException ("Invalid entry for " +
									"the steady state iterate");
						steadyMenu.setState(showSteadyState);
					} catch (IllegalArgumentException nfe) {
						JOptionPane.showMessageDialog(ModelPhasePlot.this, nfe.getMessage(),
								nfe.getClass().getSimpleName(), JOptionPane.ERROR_MESSAGE);
						showSteadyState = false;
						steadyStateIterate = 1;
					}
				} else {
					steadyStateIterate = 1;
					steadyMenu.setState(showSteadyState);
				}
				getSolution();
			}
    	});

    	nMenu.addActionListener(new java.awt.event.ActionListener () {
    		public void actionPerformed (java.awt.event.ActionEvent e) {
    			// prompt for input
    			try {
    				n = Integer.parseInt(JOptionPane.showInputDialog(ModelPhasePlot.this, 
    						"Enter the new number of iterations"));
    			} catch (NumberFormatException nfe) {
					JOptionPane.showMessageDialog(ModelPhasePlot.this, nfe.getMessage(),
							nfe.getClass().getSimpleName(), JOptionPane.ERROR_MESSAGE);
				}
    		}
    	});

    	plotVarMenu.addActionListener(new java.awt.event.ActionListener () {
    		public void actionPerformed (java.awt.event.ActionEvent e) {
    			changePlotVariables ();
    		}
    	});
    	
    	paramMenu.addActionListener(new java.awt.event.ActionListener () {
    		public void actionPerformed (java.awt.event.ActionEvent e) {
    			changeControlParameters ();
    		}
    	});
    	
    	varRangeMenu.addActionListener(new java.awt.event.ActionListener () {
    		public void actionPerformed (java.awt.event.ActionEvent e) {
    			changeVariableRange();
    		}
    	});
    	
    	paramRangeMenu.addActionListener(new java.awt.event.ActionListener () {
    		public void actionPerformed (java.awt.event.ActionEvent e) {
    			changeParameterRange();
    		}
    	});
    	
    	setParamMenu.addActionListener(new java.awt.event.ActionListener() {
    		public void actionPerformed (java.awt.event.ActionEvent e) {
    			setParameterValues();
    		}
    	});
    	
    	paramMinorMenu.addActionListener(new java.awt.event.ActionListener () {
    		public void actionPerformed (java.awt.event.ActionEvent e) {
    			changeParamMinorTick();
    		}
    	});
    	
    	paramMajorMenu.addActionListener(new java.awt.event.ActionListener () {
    		public void actionPerformed (java.awt.event.ActionEvent e) {
    			changeParamMajorTick();
    		}
    	});
    	
    	printMenu.addActionListener(new java.awt.event.ActionListener () {
    		public void actionPerformed (java.awt.event.ActionEvent e) {
    			printPhaseVariables();
    		}
    	});
    	
    	// call method to initialize phase plot
    	initializePhasePlot ();
    	
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
    	 			for (int i = steadyStateIterate; i < v.length; ++i)
        	 			g.drawLine(toPx(sx[i - 1]), toPy(sy[i - 1]), toPx(sx[i]), toPy(sy[i]));
    	 		} else {
    	 			for (int i = steadyStateIterate; i < v.length; ++i) {
        	 			g.drawLine(
        	 					toPx(v[i - 1][xIndex]), toPy(v[i - 1][yIndex]),
        	 					toPx(v[i][xIndex]), toPy(v[i][yIndex]));
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
				String xAxis = (specialPlot) ? specialLegendX : variableList[xIndex];
				String yAxis = (specialPlot) ? specialLegendY : variableList[yIndex];
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
		
        p1Label = new javax.swing.JLabel();
        p2Label = new javax.swing.JLabel();
        p2Slider = new javax.swing.JSlider();
        p1Slider = new javax.swing.JSlider();
        p1ValueLabel = new javax.swing.JLabel();
        p2ValueLabel = new javax.swing.JLabel();
        
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
					setVariableInitialCondition (toX(e.getX()), toY(e.getY()));
			}

			public void mouseReleased(MouseEvent e) {
				// TODO Auto-generated method stub
				if (e.isPopupTrigger())
					popupMenu.show(e.getComponent(), e.getX(), e.getY());
			}
        });
        
        javax.swing.GroupLayout phasePlotLayout = new javax.swing.GroupLayout(phasePlot);
        phasePlot.setLayout(phasePlotLayout);
        phasePlotLayout.setHorizontalGroup(
            phasePlotLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 498, Short.MAX_VALUE)
        );
        phasePlotLayout.setVerticalGroup(
            phasePlotLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 498, Short.MAX_VALUE)
        );

        p1Label.setFont(new java.awt.Font("Tahoma", 1, 14)); // NOI18N
        p1Label.setText(parameterList[p1Index] + ":");

        p2Label.setFont(new java.awt.Font("Tahoma", 1, 14)); // NOI18N
        p2Label.setText(parameterList[p2Index] + ":");

        p2Slider.setMajorTickSpacing(20);
        p2Slider.setMaximum(200);
        p2Slider.setMinorTickSpacing(2);
        p2Slider.setPaintTicks(true);
        p2Slider.setValue(0);
        p2Slider.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
            	evaluateP2 ();
            }
        });

        p1Slider.setMajorTickSpacing(20);
        p1Slider.setMaximum(200);
        p1Slider.setMinorTickSpacing(2);
        p1Slider.setPaintTicks(true);
        p1Slider.setValue(0);
        p1Slider.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
            	evaluateP1();
            }
        });

        p1ValueLabel.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
        p1ValueLabel.setText("0.0");

        p2ValueLabel.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
        p2ValueLabel.setText("0.0");

        javax.swing.GroupLayout jPanel2Layout = new javax.swing.GroupLayout(jPanel2);
        jPanel2.setLayout(jPanel2Layout);
        jPanel2Layout.setHorizontalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            		.addComponent(p1Label)
                    .addComponent(p2Label))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(p1Slider, javax.swing.GroupLayout.DEFAULT_SIZE, 425, Short.MAX_VALUE)
                    .addComponent(p2Slider, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            		.addComponent(p1ValueLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(p2ValueLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
            .addComponent(paramPane, 500, 500, 500)
        );
        jPanel2Layout.setVerticalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(p1Slider, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(p1Label, javax.swing.GroupLayout.PREFERRED_SIZE, 17, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(p1ValueLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 25, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(5, 5, 5)
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(p2Slider, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(p2Label)
                    .addComponent(p2ValueLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 25, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(0, 0, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGap(5, 5, 5)
                .addComponent(paramPane, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(5, 5, 5)
                )
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(10, 10, 10)
                        .addComponent(jPanel2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(phasePlot, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(phasePlot, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jPanel2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(0, 0, 0))
        );
        
        pack();
        //setLocationRelativeTo (null);
        
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
    // Variables declaration - do not modify                     
    private javax.swing.JLabel p1Label;
    private javax.swing.JLabel p2Label;
    private javax.swing.JPanel phasePlot;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JSlider p1Slider;
    private javax.swing.JSlider p2Slider;
    private javax.swing.JLabel p1ValueLabel;
    private javax.swing.JLabel p2ValueLabel;
    private javax.swing.JPopupMenu popupMenu;
    private javax.swing.JTextArea paramArea;
    
    // End of variables declaration 

}
