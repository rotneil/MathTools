package models;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.WindowEvent;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import math.tools.MatrixException;


public class LinearBifurcation extends JPanel 
{
	private static final long serialVersionUID = 1L;

	// the frame
	javax.swing.JFrame mFrame;
	private java.io.File mFile;
	
	// instance variables
	private double y0;
	private double yf;
	private double x0;
	private double xf;
	
	private int xStart = 55;
	private int yStart = 50;
	private int height = 350;
	private int width = 550;
	private int grid = 800;
	
	private int iterate = 800;		// the number of iterates to make up the vertical grid
	private int offIterate = 10000;	// the number of iterates to ignore in order to attain stability
	
	private int yDiv = 6, xDiv = 4;
	private double dt = 0.01;
	
	private java.util.ArrayList<double []> bif; 	// the linear bifurcation points
	
	// sekikaw2010 parameters
	private double e = 0.1, k1 = 0.35, k2 = k1, b0 = 0.49, b1 = b0, control;   // k3 is the control 
	private int dim = 3;				// dimension of the system;
	
	// index of the last computed bifurcation
	private int mIndex;
	
	// Constructor
	public LinearBifurcation (javax.swing.JFrame frame, double x0, double xf, double y0, double yf)
	{
		mFrame = frame;
		this.x0 = x0;
		this.xf = xf;
		this.y0 = y0;
		this.yf = yf;
		
		// instantiate the bif values
		bif = new java.util.ArrayList<double[]>();
		
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
							"Compute local Bifurcation");
					javax.swing.JMenuItem item2 = new javax.swing.JMenuItem(
							"Save Bifurcation to File");
					
					item1.addActionListener(new ActionListener() {
						@Override
						public void actionPerformed(ActionEvent e) {
							// TODO Auto-generated method stub
							computeRegionBifurcation ();
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
			while (!dialog.getFile().endsWith(".lnb")) {
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
					computeBifurcationDiagram ();
				else
					repaintPlot();
			}
			else{
				mFile.createNewFile();
				mIndex = 0;
				computeBifurcationDiagram ();
			}
		} catch (java.io.IOException e) {
			javax.swing.JOptionPane.showMessageDialog(mFrame, e.getMessage(), 
					e.getClass().getSimpleName(), javax.swing.JOptionPane.ERROR_MESSAGE);
		}
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
			double [] lnb = new double [iterate];
			for (int i = 0; i < word.length; ++i)
				lnb [i] = Double.parseDouble(word[i]);
			bif.add(lnb);
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
						for (int i = 0; i < bif.size(); ++i) {
							double [] lnb = bif.get(i);
							for (int j = 0; j < lnb.length; ++j)
								writer.write(lnb [j] + " ");
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
	private void computeRegionBifurcation ()
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
							"Enter the dimension of the region as x1,x2,y1,y2 ",
							javax.swing.JOptionPane.INFORMATION_MESSAGE);
					String [] input = range.split(",");
					
					double x1 = Double.parseDouble(input[0]);
					double x2 = Double.parseDouble(input[1]);
					double y1 = Double.parseDouble(input[2]);
					double y2 = Double.parseDouble(input[3]);
					
					// check for valid extry of x1, and x2
					if (x1 > x2 || x1 < x0 || x2 > xf)
						throw new IllegalArgumentException ("Entered horizontal range is invalid!");
					if (y1 > y2 || y1 < y0 || y2 > yf)
						throw new IllegalArgumentException ("Entered vertical range is invalid!");
					
					// launch anothe frame
					javax.swing.JFrame frame = new javax.swing.JFrame (
							"Bifurcation diagram for an autonomous system");
					frame.add(new LinearBifurcation (frame, x1, x2, y1, y2));
					frame.setSize (680, 500);
					frame.setLocationRelativeTo(mFrame);
					frame.setVisible(true);
					
					// repaint plot
					repaintPlot();
					
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
	private void computeBifurcationDiagram () {
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {
			public void run () {
				// compute the arrayList
				
				// iterate through the grid-- assigning the value of the control
				// parameter and computing the bifurcation of the stable points
				for ( ; mIndex < grid + 1; ++mIndex) {
					// k3 is the control paramter for sekikawa2010
					control = x0 + mIndex * (xf - x0) / grid;
					
					bif.add(computeBifurcation ());
					
					// show the plotted shape
					repaintPlot();
				}
				repaintPlot();
				saveToFile();
			}
		});
		executor.shutdown();
		repaint();
	}
	
	// method to compute the specific bifurcation diagram for a given value of the
	// control parameter
	private double [] computeBifurcation ()
	{	
		// set the system to an initial condition
		double [] v = new double [dim];
		
		for (int i = 0; i < dim; ++i)
			v[i] = random (y0, yf);
		
		// perform trade-off iteration to attain stability
		for (int i = 0; i < offIterate; ++i)
			v = rungeKutta (v);
		
		return getMinima (v);
	}
	
	// method to solve equation
	private double [] getMinima (double [] v0)
	{
		// local variables
		double minima [] = new double [iterate];
		minima[0] = control;
		
		double v1 [] = new double [v0.length];
		double v2 [] = new double [v0.length];
		
		for (int i = 1; i < minima.length; ) {
			// evaluate new coordinate
			v1 = rungeKutta (v0);
			v2 = rungeKutta (v1);
			
			// assume that the minimum point is at v1
			if (v1[0] < v0[0] && v1[0] < v2[0])		// as applicable to sekikawa2010
			{
				// assign maximum value
				minima [i] = v1[0];
				
				// increment counter
				++i;
			}
			// set the new initial value
			v0 = v1;
		}
		
		return minima;
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
	
	// equation definition
	private double [] f (double [] v)
	{
		int l = v.length;
		double [] vv = new double [l];
		
		// sekikawa2010
		vv[0] = (v[0] * (1 - sq (v[0])) + v[1] + v[2]) / e;
		vv[1] = -v[0] - k1 * v[1] + b0;
		vv[2] = control * (-v[0] - k2 * v[2] + b1);
		
		return vv;
	}
	
	
	// UTILITY METHODS
	private double sq (double x) {return x * x; }
	private double random (double x1, double x2) {return x1 + (x2 - x1) * Math.random(); }
	
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
	
	// method to paint on the panel
	public void paintComponent (Graphics g)
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
		g.drawString("k3", (int)(xStart + 0.5 * width), yStart + height + 35);
		g.drawString("x", 15, (int)(yStart + 0.5 * height));
		
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
		
		
		// draw plot on the panel
		int px, py;
		for (int i = 0; i < bif.size(); ++i) {
			// get the linear bifurcation
			double lnb [] = bif.get(i);
			px = toPx(lnb [0]);
			
			// iterate through the vertical points
			for (int j = 1; j < lnb.length; ++j) {
				py = yStart + (int) ((yf - lnb [j]) / (yf - y0) * height);
				g.drawOval (px, py, 1, 1);
			}
		}		
	}
	
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
						"Bifurcation Diagram for an autonomous system");
				frame.add(new LinearBifurcation (frame, 0.35, 0.75, -0.85, 0.45));
				frame.setSize (680, 500);
				//frame.setLocationRelativeTo(null);
				frame.setVisible(true);
			}
		});
	}
}
