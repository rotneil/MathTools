package models;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.WindowEvent;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import math.tools.Matrix;
import math.tools.MatrixException;


public class LyapunovExponents extends javax.swing.JPanel
{
	private static final long serialVersionUID = 1L;
	
	// instance variables
	private javax.swing.JFrame mFrame;
	
	private int xStart = 55;
	private int yStart = 50;
	private int height = 350;
	private int width = 550;
	
	// plot dimension for sekikawa2010
	/*
	private double x0 = 0.54;	// magnified
	private double xf = 0.57;	// magnified
	private double y0 = -0.4;	// magnified
	private double yf = 0.1;	// magnified
	*/
	private double x0 = 0.35;
	private double xf = 0.75;
	private double y0 = -1.2;
	private double yf = 0.2;
	
	private int xDiv = 4;
	//private int yDiv = 5;	// magnified
	private int yDiv = 7;
	
	/**
	 ** plot dimension for ueta2007
	private double x0 = 1.0;
	private double xf = 1.5;
	private double y0 = -0.8;
	private double yf = 0.2;
	
	private int xDiv = 5;
	private int yDiv = 5;
	 */
	
	
	/**
	** plot dimension for Ueta2003
	**
	private double x0 = 1.8;
	private double xf = 2.15;
	private double y0 = -0.5;
	private double yf = 0.1;
	
	private int xDiv = 7;
	private int yDiv = 6;
	*/
	
	private int grid = 1000;
	private double dt = 0.01;
	
	private java.util.ArrayList<double[]> lyap;
	
	// term definition for sekikaw2010
	private double e = 0.1, k1 = 0.35, k2 = k1, b0 = 0.49, b1 = b0;	// k3 is the control param
	private int offIterate = 10000;
	
	// term definition for ueta2004
	//private double g = 1.637, k2 = 0.56, d = 0.337; // k1 is the control param
	// private int offIterate = 500000;
	
	// term definition for ueta2003
	//private double g = 1.6365, k1 = 0.653, k2 = 1.05, d;	// d is the control param
	// private int offIterate = 100000;
	
	private int dim;
	private int mIndex;
	
	// the files into which to save updates
	private java.io.File mFile;
	
	// constructor
	public LyapunovExponents (javax.swing.JFrame frame)
	{
		// instantiate the lyapunov exponents
		mFrame = frame;
		dim = 3;
		lyap = new java.util.ArrayList<>();
		
		// set the window closing event for frame
		mFrame.addWindowListener(new java.awt.event.WindowAdapter () {
			@Override
			public void windowOpened(WindowEvent e) {
				// TODO Auto-generated method stub
				ExecutorService executor = Executors.newCachedThreadPool();
				executor.execute(new Runnable () {
					@Override
					public void run () {
						launchFileDialog();
					}
				});
				executor.shutdown();
			}
			
			@Override
			public void windowClosed(WindowEvent e) {
				saveToFile();
			}
		});
		
		mFrame.addMouseListener(new java.awt.event.MouseAdapter () {
			@Override
			public void mouseReleased(MouseEvent e) {
				// TODO Auto-generated method stub
				if (e.isPopupTrigger()) {
					javax.swing.JPopupMenu popup = new javax.swing.JPopupMenu();
					javax.swing.JMenuItem item1 = new javax.swing.JMenuItem (
							"Compute extra exponents");
					javax.swing.JMenuItem item2 = new javax.swing.JMenuItem("Save Exponents to File");
					
					item1.addActionListener(new ActionListener() {
						@Override
						public void actionPerformed(ActionEvent e) {
							// TODO Auto-generated method stub
							computeExtraExponent();
						}
					});
					
					item2.addActionListener(new ActionListener() {
						@Override
						public void actionPerformed(ActionEvent e) {
							// TODO Auto-generated method stub
							saveToFile();
						}
					});
					popup.add(item1);
					popup.add(item2);
					popup.show(e.getComponent(), e.getX(), e.getY());
				}
			}
		});
	}
	
	/**
	 * This method launches the file dialog
	 */
	private void launchFileDialog ()
	{		
		// prompt the user to save to file or load from it
		java.awt.FileDialog dialog = new java.awt.FileDialog (mFrame,
				"Load exponents from the selected file or save to it");
		dialog.setMultipleMode(false);
		dialog.setVisible(true);
		
		// check that the dialog was not cancelled
		if (dialog.getFile() != null) {
			while (!dialog.getFile().endsWith(".lxp")) {
				javax.swing.JOptionPane.showMessageDialog(mFrame, "Wrong file type selected",
						"WRONG SELECTION", javax.swing.JOptionPane.ERROR_MESSAGE);
				dialog.setVisible(true);
			}
			if (dialog.getFile() != null)
				createFile(dialog.getDirectory(), dialog.getFile());
			else
				System.exit(0);
		}
		else
			// shutdown the application
			System.exit(0);
	}
	
	/**
	 * This method creates the file to be used for storing exponents for
	 * future retrieving of the exponents
	 * @param directory
	 * @param file
	 */
	private void createFile (String directory, String file)
	{
		// create file
		try {
			mFile = new java.io.File(directory, file);
			
			if (mFile.exists()) {
				mIndex = loadFromFile();
				if (mIndex < grid)
					computeExponents ();
				else
					repaintPlot();
			}
			else{
				mFile.createNewFile();
				mIndex = 0;
				computeExponents ();
			}
		} catch (java.io.IOException e) {
			javax.swing.JOptionPane.showMessageDialog(mFrame, e.getMessage(), 
					e.getClass().getSimpleName(), javax.swing.JOptionPane.ERROR_MESSAGE);
		}
	}
	
	/**
	 * This method reads the content of the specified file into a 
	 * two-dimensional double-precision array lyap.
	 * @return
	 * @throws java.io.IOException
	 */
	private int loadFromFile () throws java.io.IOException
	{
		if (mFile == null)
			throw new java.io.IOException ("NullPointerException. The chosen file is empty.");
		
		java.io.BufferedReader reader = new java.io.BufferedReader (
				new java.io.FileReader(mFile));
		int index = 0;
		while (reader.ready()) {
			String [] word = reader.readLine().split(" ");
			double [] lxp = new double [dim + 1];
			for (int i = 0; i < word.length; ++i)
				lxp[i] = Double.parseDouble(word[i]);
			lyap.add(lxp);
			++index;
		}
		reader.close();
		
		return index;
	}
	
	/**
	 * This method writes the content of the two-dimensional array into a file starting
	 * from index zero to the specified index of the array.
	 * @throws java.io.IOException
	 */
	private void saveToFile ()
	{
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {
			@Override
			public void run () {
				try {
					// check if mWriter is not null
					if (mFile != null) {
						// instantiate a random access file writer
						java.io.FileWriter writer = new java.io.FileWriter(mFile);
						for (int i = 0; i < lyap.size(); ++i) {
							double [] lxp = lyap.get(i);
							for (int j = 0; j < lxp.length; ++j)
								writer.write(lxp [j] + " ");
							writer.write("\n");
						}
						
						writer.flush();
						writer.close();
					}
				} catch (java.io.IOException e) {
					javax.swing.JOptionPane.showMessageDialog(mFrame, e.getMessage(), 
							e.getClass().getSimpleName(), javax.swing.JOptionPane.ERROR_MESSAGE);
				}
			}
		});
	}
	
	// method to start the computation of lyapunov exponents for extra features
	private void computeExtraExponent ()
	{
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {
			@Override
			public void run() {
				try {
					// check that xFile is not null
					if (mFile == null)
						launchFileDialog();
					
					// prompt the user to choose the range of the control param to compute
					String range = javax.swing.JOptionPane.showInputDialog(mFrame, 
							"Enter the range as x1,x2 of the control parameter",
							"Extra Exponent Computation Dialog",
							javax.swing.JOptionPane.INFORMATION_MESSAGE);
					String [] input = range.split(",");
					
					double x1 = Double.parseDouble(input[0]);
					double x2 = Double.parseDouble(input[1]);
					
					// check for valid extry of x1, and x2
					if (x1 > x2 || x1 < x0 || x2 > xf)
						throw new IllegalArgumentException ("Entered range is invalid");
					
					double dx = (xf - x0) / grid;
					
					while (x1 <= x2) {
						// add to the arraylist of exponent extra
						lyap.add(computeLyapunov(x1));
						
						// repaint
						repaintPlot();
						x1 += dx;
					}
					saveToFile();
					
				} catch (MatrixException e) {
					javax.swing.JOptionPane.showMessageDialog(mFrame, e.getMessage(),
							"Inverse Computation Error",
							javax.swing.JOptionPane.ERROR_MESSAGE);
				} catch (IllegalArgumentException e) {
					javax.swing.JOptionPane.showMessageDialog(mFrame, e.getMessage(),
							"Extra Exponent Computation Error",
							javax.swing.JOptionPane.ERROR_MESSAGE);
				}
			}
		});
		executor.shutdown();
	}
	
	// method to start the computation of lyapunov exponents
	private void computeExponents () {
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {
			public void run () {
				// compute the exponents
				try {
					// control parameter
					double control;
					
					// iterate through the grid-- assigning the value of the control
					// parameter and computing the lyapunov exponents
					for ( ; mIndex < grid; ++mIndex) {
						// k1 is the control parameter for the Ueta2004
						control = x0 + mIndex * (xf - x0) / grid;
						
						// d is the control parameter for Ueta2003
						// d = x0 + i * (xf - x0) / grid;
						lyap.add(computeLyapunov(control));
						
						// show the plotted shape
						repaintPlot();
					}
					repaintPlot();
					saveToFile();
				} catch (MatrixException e) {
					javax.swing.JOptionPane.showMessageDialog(mFrame, e.getMessage(), 
							e.getClass().getSimpleName(), javax.swing.JOptionPane.ERROR_MESSAGE);
				}
			}
		});
		executor.shutdown();
		repaint();
	}
	
	// method to repaint plot
	private void repaintPlot ()
	{
		javax.swing.SwingUtilities.invokeLater(new Runnable (){
			public void run () {
				repaint();
			}
		});
	}
	
	/**
	 * This method computes the lyapunov for the current iterate
	 * @param 
	 * @return
	 */
	private double [] computeLyapunov (double k1) throws MatrixException
	{
		// local variables
		int n = 150000;
		double [] h = new double [dim + 1];
		h[0] = k1;			// save the current control parameter value
		
		// set initial variation
		double [] v = setInitialVariation();
		
		//Perform i-iterate of the orbit in order to reach the attractor
		for (int i = 0; i < offIterate; ++i)
			v = rungeKutta (v, k1);
		
		
		double [] [] w = Matrix.getIdentity (dim);
		double [] [] z, Df;
		
		// perform Lyapunov computation via Gram-Schmidt Orthognalization
		for (int i = 0; i < n; ++i) {
			// reset the variational variables
			resetVariationVariables (v);
			
			// Iterate the orbit in order to make a Time-T map
			for (int j = 0; j < 100; ++j)
				v = rungeKutta (v, k1);
			
			// get the ellipsoid
			Df = extractEllipsoid (v);
			
			z = Matrix.product (Df, w);
			
			// orthogonalize z
			w = Matrix.getOrthogonalVector (z);
			
			// get the expansion and convert to lyapunov exponents
			double [] len = Matrix.getModulus (w);
			for (int j = 0; j < len.length; ++j)
				h [j + 1] += Math.log(len[j]);
			
			// normalize the vector
			Matrix.normalize (w);
		}
		
		for (int t = 1; t < h.length; ++t)
			h[t] /= n;
		
		return h;
	}
	
	// this method resets the variation matrix for futher iterates
	private void resetVariationVariables (double [] v)
	{
		// set all variables from dim to the end of array to zero
		for (int i = dim; i < v.length; ++i)
			v[i] = 0.0;
		for (int i = 0; i < dim; ++i)
			v[(i + 1) * dim + i] = 1.0;
	}
	
	// the method to extract the ellipsoid from the matrix v
	private double [] [] extractEllipsoid (double [] v)
	{
		// local variables
		double [] [] ellipsoid = new double [dim][dim];
		/*
		new double [] [] {
				{v[4], v[8],  v[12], v[16]}, 
				{v[5], v[9],  v[13], v[17]},
				{v[6], v[10], v[14], v[18]},
				{v[7], v[11], v[15], v[19]}
		};*/
		// start the extraction
		for (int i = 0; i < dim; ++i)
			for (int j = 0; j < dim; ++j)
				ellipsoid[i][j] = v[i + (j + 1) * dim];
		
		return ellipsoid;
	}
	
	// the method that sets the initial variation
	private double [] setInitialVariation ()
	{
		// local variable
		double [] init = new double [(1 + dim) * dim];
		
		for (int i = 0; i < dim; ++i)
			init [i] = random ();
		for (int i = 0; i < dim; ++i)
			init [(i + 1) * dim + i] = 1.0;
		
		return init;
	}
	
	/**
	 * The function of the time continuous differential equation
	 * @param v
	 * @return
	 */
	private double [] f (double [] v, double control)
	{
		// declare local variable
		double [] vv = new double [v.length];
		
		// expressions for sekikawa2010
		vv[0] = (v[0] * (1 - sq(v[0])) + v[1] + v[2]) / e;
		vv[1] = -v[0] - k1 * v[1] + b0;
		vv[2] = control * (-v[0] - k2 * v[2] + b1);
		
		for (int i = 0; i < 9; i += dim) {
			vv[i + 3] = (1 - 3 * sq(v[0])) / e * v[i + 3] + v[i + 4] / e + v[i + 5] / e;
			vv[i + 4] = -v[i + 3] - k1 * v[i + 4];
			vv[i + 5] = -control * v[i + 3] - k2 * control * v[i + 5];
		}
		
		/**
		 ** The expression below represents the function for ueta2004
		 **
		// the definition as relevant to Ueta2004
		vv[0] = -v[1] + Math.tanh(g * v[0]) - d * (v[0] - v[2]);
		vv[1] = v[0] - control * v[1];
		vv[2] = -v[3] + Math.tanh(g * v[2]) - d * (v[2] - v[0]);
		vv[3] = v[2] - k2 * v[3];
		
		for (int i = 0; i < 16; i += dim) {
			vv[i + 4] = (g / sq(Math.cosh(g*v[0])) - d) * v[i + 4] - v[i + 5] + d * v[i + 6];
			vv[i + 5] = v[i + 4] - control * v[i + 5];
			vv[i + 6] = d * v[i + 4] + (g / sq(Math.cosh(g * v[2])) - d) * v[i + 6] - v[i + 7];
			vv[i + 7] = v[i + 6] - k2 * v[i + 7];
		}
		*/
		
		/**
		 ** The expressions below represent the function for Ueta2003
		 **
		double n = 1.0 / (1 + d * k1);
		
		vv[0] = -v[1] + Math.tanh(g*v[0]);
		vv[1] = v[0] - control * v[1] + d * k1 * n * (control * v[1] - v[2]);
		vv[2] = -v[3] + Math.tanh(g * v[2]) + d * n * (control * v[1] - v[2]);
		vv[3] = v[2] - k2 * v[3];
		
		for (int i = 0; i < 16; i += 4) {
			vv[i + 4] = g / sq(Math.cosh(g * v[0])) * v[i + 4] - v[i + 5];
			vv[i + 5] = v[i + 4] + (d * sq(k1) * n - control) * v[i + 5] - 
				d * control * n * v[i + 6];
			vv[i + 6] = control * d * n * v[i + 5] + 
					(g / sq(Math.cosh(g * v[2])) - d * n) * v[i + 6] - v[i + 7];
			vv[i + 7] = v[i + 6] - k2 * v[i + 7];
		}
		*/
		return vv;
	}
		
    // Runge Kutta fourth order numerical integeration
 	private double [] rungeKutta (double [] v, double k1)
 	{
 		int l = v.length;
 		double [] c1, c2, c3, c4;
 		
 		// initialize the intermediate steps
 		c1 = new double [l];
 		c2 = new double [l];
 		c3 = new double [l];
 		c4 = new double [l];
 		
 		c1 = f(v, k1);
 		
 		for (int i = 0; i < l; ++i)
 			c2[i] = v[i] + dt * c1[i] / 2;
 		c2 = f(c2, k1);
 		
 		for (int i = 0; i < l; ++i)
 			c3[i] = v[i] + dt * c2[i] / 2;
 		c3 = f(c3, k1);
 		
 		for (int i = 0; i < l; ++i)
 			c4[i] = v[i] + dt * c3[i];
 		c4 = f(c4, k1);
 		
 		for (int i = 0; i < l; ++i)
 			c1[i] = v[i] + dt * (c1[i] + 2 * (c2[i] + c3[i]) + c4[i]) / 6.0;
 		
 		return c1;
 	}
 	
	
	@Override
	public void paintComponent (java.awt.Graphics g)
	{
		super.paintComponent(g);
		
		// turn the background to white
		setBackground(java.awt.Color.white);
		
		// show the draw area
		g.drawRect(xStart, yStart, width, height);
		
		// show the legend
		g.drawString("" + yf, xStart - 25, yStart + 5);
		g.drawString("" + y0, xStart - 25, yStart + height + 5);
		g.drawString("" + x0, xStart - 5, yStart + height + 15);
		g.drawString("" + xf, xStart + width - 10, yStart + height + 15);
		
		// show the plot variables
		g.drawString("d", (int)(xStart + 0.5 * width), yStart + height + 25);
		g.drawString("v", 15, (int)(yStart + 0.5 * height));
		
		// show major marks for the y-axis
		int hh;
		for (int i = 1; i < yDiv; ++i) {
			hh = yStart + (int) (i * height / yDiv) + 1;
			g.drawLine(xStart, hh, xStart + 5, hh);
			g.drawLine(xStart + width - 5, hh, xStart + width, hh);
			g.drawString("" + todp ((yf - i * (yf - y0) / yDiv), 1), xStart - 25, hh + 5);
		}
		
		// show major marks for the x-axis
		for (int i = 1; i < xDiv; ++i) {
			hh = xStart + (int) (i * width / xDiv);
			g.drawLine(hh, yStart, hh, yStart + 5);
			g.drawLine(hh, yStart + height - 5, hh, yStart + height);
			g.drawString("" + todp((x0 + i * (xf - x0) / xDiv), 2), hh - 10, yStart + height + 15);
		}
		
		
		// plot the scatter plot of the exponents
		for (int i = 0; i < lyap.size() - 1; ++i) {
			double [] lxp1 = lyap.get(i);
			double [] lxp2 = lyap.get(i + 1);
			int px1 = toPx(lxp1[0]);
			int px2 = toPx(lxp2[0]);
			
			for (int j = 1; j < lxp1.length; ++j) {
				int py1 = yStart + (int) ((yf - lxp1[j]) / (yf - y0) * height);
				int py2 = yStart + (int) ((yf - lxp2[j]) / (yf - y0) * height);
				g.drawLine (px1, py1, px2, py2);
			}
		}
		
	}   // end method paintComponent
	
	/**
	 * This method converts a value x to pixels
	 * @param x
	 * @return
	 */
	private int toPx(double x)
	{
		return (int) (xStart + (x - x0) / (xf - x0) * width);
	}
	
	/**
	 * Mathod to compute the square of a number x
	 * @param x
	 * @return
	 */
	private double sq (double x) {return x * x; }
	private double random () {return -2.0 + 4.0 * Math.random(); }
	
	
	/**
	 * This method converts a double-value to a certain decimal place
	 * @param value
	 * @param dp
	 * @return
	 */
	private double todp (double value, int dp)
	{
		double pow = Math.pow(10, dp);
		int v = (int) (value * pow + Math.signum(value) * 0.5);
		return 1.0 * v / pow;
	}
	
	/**
	 * This method execute the LyapunovExponent Application
	 * @param args
	 */
	public static void main (String [] args)
	{
		javax.swing.SwingUtilities.invokeLater(new Runnable () {
			@Override
			public void run () {
				// instantiate the frame upon which to display the exponents
				javax.swing.JFrame frame = new javax.swing.JFrame (
						"Lyapunove Exponents for an autonomous system");
				frame.add(new LyapunovExponents(frame));
				frame.setSize (680, 500);
				//frame.setLocationRelativeTo(null);
				frame.setVisible(true);
			}
		});
	}
}
