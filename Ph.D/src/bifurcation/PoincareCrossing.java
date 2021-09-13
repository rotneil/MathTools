package bifurcation;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.JFrame;

public class PoincareCrossing extends javax.swing.JPanel
{
	// instance variables
	private int xStart = 55;
	private int yStart = 20;
	private int height = 450;
	private int width = 450;
	
	// plot dimension
	private double x0;
	private double xf;
	private double y0;
	private double yf;
	private int grid = 1000;
	
	private double dt = 0.01;
	private double [] x, y;
	private static final int n = 10000000;
	private int counter = 0;
	private boolean backward = false;
	private java.util.ArrayList <Point> fixedPoint;
	
	private double g = 0.0;
	private double k = 0.0;
	private double dp = 0.0001;
	private javax.swing.JFrame mFrame;
	
	// constructor
	public PoincareCrossing (javax.swing.JFrame f, 
			final double x0, final double xf, final double y0, final double yf)
	{
		// instantiate the instance data
		this.x0 = x0;
		this.xf = xf;
		this.y0 = y0;
		this.yf = yf;
		this.mFrame = f;
		
		x = new double [n];
		y = new double [n];
		fixedPoint = getFixedPoint (g, k);
		solveTrajectory();
		
		// set plot background to white
		setBackground (java.awt.Color.white);
		
		// listen to mouse clicks
		addMouseListener (new MouseAdapter () {
        	public void mouseClicked (MouseEvent e) {
        		// retrieve the points
        		x[0] = toX (e.getX());
        		y[0] = toY (e.getY());
        		counter = 0;
        		solveTrajectory();
        	}
        });
		
		// listen to key event
		mFrame.addKeyListener(new java.awt.event.KeyAdapter () {
			public void keyTyped(KeyEvent e) {
				// TODO Auto-generated method stub
				switch (e.getKeyChar()) {
					case 'z': case 'Z':
						
						break;
					case 'i': case 'I':
						// prompt for initial condition
						try {
							x[0] = Double.parseDouble(javax.swing.JOptionPane.showInputDialog(
									mFrame, "Enter initial value for x", "Input Dialog",
									javax.swing.JOptionPane.INFORMATION_MESSAGE));
							y[0] = Double.parseDouble(javax.swing.JOptionPane.showInputDialog(
									mFrame, "Enter initial value for y", "Input Dialog",
									javax.swing.JOptionPane.INFORMATION_MESSAGE));
							solveTrajectory();
							counter = 0;
						} catch (Exception ex) {
							javax.swing.JOptionPane.showMessageDialog(mFrame, 
									ex.getClass().getSimpleName(), ex.getMessage(), 
									javax.swing.JOptionPane.INFORMATION_MESSAGE);
						}
						break; 
					case 'x' : case 'X':
						try {
							g = Double.parseDouble(javax.swing.JOptionPane.showInputDialog(
									mFrame, "Enter the value of g: ", "Input Dialog",
									javax.swing.JOptionPane.INFORMATION_MESSAGE));
							k = Double.parseDouble(javax.swing.JOptionPane.showInputDialog(
									mFrame, "Enter the value of k: ", "Input Dialog",
									javax.swing.JOptionPane.INFORMATION_MESSAGE));
							fixedPoint = getFixedPoint (g, k);
							solveTrajectory();
							counter = 0;
						} catch (Exception ex) {
							javax.swing.JOptionPane.showMessageDialog(mFrame, 
									ex.getClass().getSimpleName(), ex.getMessage(), 
									javax.swing.JOptionPane.INFORMATION_MESSAGE);
						}
						break;
					case 'c' : case 'C':
						if (!backward && counter < n)
							++counter;
						else if (backward && counter > 0)
							--counter;
						repaint();
						break;
					case 'v' : case 'V':
						if (!backward && counter < n - 10)
							counter += 10;
						else if (backward && counter > 10)
							counter -= 10;
						repaint();
						break;
					case 'b' : case 'B':
						if (!backward && counter < n - 100)
							counter += 100;
						else if (backward && counter > 100)
							counter -= 100;
						repaint();
						break;
					case 'n' : case 'N':
						if (!backward && counter < n - 1000)
							counter += 1000;
						else if (backward && counter > 1000)
							counter -= 1000;
						repaint();
						break;
					case 'm' : case 'M':
						if (!backward && counter < n - 10000)
							counter += 10000;
						else if (backward && counter > 10000)
							counter -= 10000;
						repaint();
						break;
					case '-':
						backward = true;
						if (counter > 0)
							--counter;
						repaint();
						break;
					case '+': case '=':
						backward = false;
						if (counter < n)
							++counter;
						repaint();
						break;
					case 'g' : case 'G':
						if (!backward && g < 2.0)
							g += dp;
						else if (backward && g > 0)
							g -= dp;
						fixedPoint = getFixedPoint (g, k);
						solveTrajectory();
						counter = 0;
						repaint();
						break;
					case 'k' : case 'K':
						if (!backward && k < 2.0)
							k += dp;
						else if (backward && k > 0)
							k -= dp;
						fixedPoint = getFixedPoint (g, k);
						solveTrajectory();
						counter = 0;
						repaint();
						break;
					case 'p' : case 'P':
						try {
							dp = Double.parseDouble(
									javax.swing.JOptionPane.showInputDialog(mFrame,
										"Enter the incremental value dp:",
										"Input Dialog",
										javax.swing.JOptionPane.INFORMATION_MESSAGE));
						} catch (Exception ex) {
							javax.swing.JOptionPane.showMessageDialog(mFrame, 
									ex.getClass().getSimpleName(), ex.getMessage(),
									javax.swing.JOptionPane.ERROR_MESSAGE);
						}
						break;
					default:
				}
			}
		});
	}
	
	// method to solve the nonlinear equation according to the supplied variables
	private void solveTrajectory ()
	{
		// local variable
		double [] vv = new double [2];
		
		// iterate orbit n-times
		for (int i = 1; i < x.length; ++i) {
			vv = rungeKutta (g, k, new double [] {x [i - 1], y [i - 1]});
			x[i] = vv[0];
			y[i] = vv[1];
		}
		repaint();
	}
	
	// method to tell whether there is a crossing or not
	private boolean getCrossing ()
	{
		if (counter <= 0)
			return false;
		
		// local variables
		double section = getPoincareSection ();
		double [] v, vv;
		
		// set the values of v0 and v1
		v = rungeKutta (g, k, new double [] {x[counter - 1], y[counter - 1]});
		vv = rungeKutta (g, k, new double [] {x[counter], y[counter]});
				
		if (v[1] < v[0] * Math.tan(section) && vv[1] >= vv[0] * Math.tan(section))
			return true;
		
		return false;
	}
	
	// method to show the difference between Poincare section and trajectory
	private double getSectionDifference ()
	{
		// local variable
		double section = getPoincareSection();
		
		// declare the y-value of the orbit on Poincare section
		double py = x[counter] * Math.tan(section);
		
		// the difference
		return (py - y[counter]);
	}
	
	// method to return Poincare section angle
	private double getPoincareSection ()
	{

		// set the section angle
		if (g * k <= 1.0)
			return Math.PI / 4.0;
		else {
			// get the fixed points
			java.util.ArrayList<Point> fixedPoint = getFixedPoint (g, k);
			Point p1 = fixedPoint.get(1);
			Point p2 = fixedPoint.get(2);
			
			// get the coordinates for evaluation of tan
			double x1, x2, y1, y2;
			x1 = p1.getX();
			x2 = p2.getX();
			y1 = p1.getY();
			y2 = p2.getY();
			
			return Math.atan((y2 - y1) / (x2 - x1));
		}
	}
	
    // Runge Kutta fourth order numerical integeration
 	private double [] rungeKutta (double g, double k, double [] v)
 	{
 		int l = v.length;
 		double [] c1, c2, c3, c4;
 		
 		// initialize the intermediate steps
 		c1 = new double [l];
 		c2 = new double [l];
 		c3 = new double [l];
 		c4 = new double [l];
 		
 		c1 = f(g, k, v);
 		
 		for (int i = 0; i < l; ++i)
 			c2[i] = v[i] + dt * c1[i] / 2;
 		c2 = f(g, k, c2);
 		
 		for (int i = 0; i < l; ++i)
 			c3[i] = v[i] + dt * c2[i] / 2;
 		c3 = f(g, k, c3);
 		
 		for (int i = 0; i < l; ++i)
 			c4[i] = v[i] + dt * c3[i];
 		c4 = f(g, k, c4);
 		
 		for (int i = 0; i < l; ++i)
 			c1[i] = v[i] + dt * (c1[i] + 2 * (c2[i] + c3[i]) + c4[i]) / 6.0;
 		
 		return c1;
 	}
 	
 	// equation definition
 	private double [] f (double g, double k, double [] v)
 	{
 		int l = v.length;
 		double [] vv = new double [l];
 		
 		// the Ueta BVP system
 		vv[0] = -v[1] + Math.tanh(g * v[0]);
 		vv[1] = v[0] - k * v[1];
 		
 		return vv;
 	}
	
	/**
	 * This method computes the fixed points of the autonomous system
	 * @param g
	 * @param k
	 * @return
	 */
    private java.util.ArrayList<Point> getFixedPoint (double g, double k)
    {
    	// declare local variables
    	java.util.ArrayList<Point> fixedPoint = new java.util.ArrayList<Point>();
    	
    	// the origin is a constant fixed point
        fixedPoint.add(new Point (0.0, 0.0));
        double m = 0.0;
        
    	// the other points
        if (g * k > 1) {
        	m = Math.sqrt(3.0 * (g * k - 1.0) / (k * g * g * g));
        	
	        // use Newton-Raphson to refine the equilibrium point
	 		m = refinePoint (g, k, m);
	 		
	 		double n = m / k;
	 		fixedPoint.add(new Point (-m, -n));
	 		fixedPoint.add(new Point (m, n));
        }
        
        return fixedPoint;
    }
    
    /**
     * This method refines the rough estimate of the fixed point
     * @param g
     * @param k
     * @param x0
     * @return
     */
    private double refinePoint (double g, double k, double x0)
 	{
 		// local variable
 		double fx = Math.tanh(g * x0) - x0 / k;
 		double ffx = g * Math.pow(Math.cosh(g * x0), -2.0) - 1.0 / k;
 		double x1 = x0 - fx / ffx;
 		
 		while (Math.abs(x1 - x0) >= 0.0001) {
 			// the new iterate
 			x0 = x1;
 			fx = Math.tanh(g * x0) - x0 / k;
 			ffx = g * Math.pow(Math.cosh(g * x0), -2.0) - 1.0 / k;
 			x1 = x0 - fx / ffx;
 		}
 		
 		return x1;
 	}
 	
 	// PAINT COMPONENT
 	public void paintComponent (Graphics g)
 	{
		super.paintComponent(g);
	
 		// draw the inner rectangle
 		g.drawRect(xStart, yStart, width, height);

 		// draw a straight line diagonally or across the fixed points
 		if (fixedPoint.size() > 1) {
 			Point p1, p2;
 			double x1, x2, y1, y2;
 			p1 = fixedPoint.get(1);
 			p2 = fixedPoint.get(2);
 			x1 = p1.getX(); y1 = p1.getY();
 			x2 = p2.getX(); y2 = p2.getY();
 			double slope = Math.atan((y2 - y1) / (x2 - x1));
 			int yy1 = toPy (y1 - (x1 - x0) * Math.tan(slope));
 			int yy2 = toPy (y1 + (xf - x1) * Math.tan(slope));
 			g.drawLine(xStart, yy1, xStart + width, yy2);
 		} 
 		else
 			g.drawLine(xStart, yStart + height, xStart + width, yStart);
 		
 		// mark the equilibrium points
 		g.setColor(Color.red);
 		
 		// origin
 		for (Point p : fixedPoint)
 			drawCross (g, p.getX(), p.getY());
	
 		// plot the phase state
 		g.setColor(Color.black);
 		for (int i = 1; i <= counter; ++i) {
 			g.drawLine(
				(int)(xStart + (x[i - 1] - x0) / (xf - x0) * width),
				(int)(yStart + (y[i - 1] - yf) / (y0 - yf) * height),
				(int)(xStart + (x[i] - x0) / (xf - x0) * width),
				(int)(yStart + (y[i] - yf) / (y0 - yf) * height));
 		}
 		
 		// show the time signature and state
 		g.drawString("Time: " + counter, xStart, 490);
 		g.drawString("x = " + todp(x[counter], 8), (int) (xStart + width * 0.25), 490);
 		g.drawString("y = " + todp(y[counter], 8), (int) (xStart + width * 0.70), 490);
 		
 		// perform Poincare-crossing analysis
 		g.drawString("Crossing: " + (getCrossing () ? "Yes" : "No"), xStart, 510);
 		g.drawString("dP = " + todp(getSectionDifference (), 8), (int)(xStart + width * 0.25), 510);
 		g.drawString("g = " + todp(this.g, 4), (int) (xStart + width * 0.70), 510);
 		g.drawString("k = " + todp(k, 4), (int) (xStart + width * 0.85), 510);
 	}
 	
 	// method to draw a cross
 	private void drawCross (Graphics g, double x, double y) 
 	{
 		// draw the horizontal line
 		g.drawLine(
			(int)(xStart + (x - 0.1 - x0) / (xf - x0) * width),
			(int)(yStart + (y - yf) / (y0 - yf) * height),
			(int)(xStart + (x + 0.1 - x0) / (xf - x0) * width),
			(int)(yStart + (y- yf) / (y0 - yf) * height));
		// draw the vertical line
		g.drawLine(
			(int)(xStart + (x - x0) / (xf - x0) * width),
			(int)(yStart + (y - 0.1 - yf) / (y0 - yf) * height),
			(int)(xStart + (x - x0) / (xf - x0) * width),
			(int)(yStart + (y + 0.1 - yf) / (y0 - yf) * height));
 		
 	}
	
	/**
	 * This method repaints components
	 */
	public void updatePaint ()
	{
		javax.swing.SwingUtilities.invokeLater(new Runnable () {
			public void run () {
				repaint ();
			}
		});
	}
	
	private double sq (double x) {return x * x;}
	
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

	private double todp (double d, int p)
	{
		int v = (int) (d * Math.pow(10, p) + 0.5);
		return (v / Math.pow(10, p));
	}
	// MY POINT DEFINITION
	private class Point
	{
		// instance variables
		double x = 0.0;
		double y = 0.0;
		java.text.DecimalFormat df = new java.text.DecimalFormat(" 0.0000;-0.0000");
		
		// constructor 
		public Point () { new Point (0.0, 0.0); }
		
		public Point (double x, double y)
		{
			setX(x);
			setY(y);
		}
		
		// SET METHODS
		public void setX (double x) {this.x = x;}
		public void setY (double y) {this.y = y;}
		
		// GET METHODS
		public double getX () { return this.x;}
		public double getY () { return this.y;}
		
		// method to string
		public String toString () {
			return "[" + df.format(this.x) + ", " + df.format(this.y) + "]";
		}
	}

	// CLASS COMPLEX NUMBER
	private class ComplexNumber
	{
		// instance variables
		private double r;
		private double im;
		
		// argument constructor
		public ComplexNumber (double real, double im)
		{
			this.r = real;
			this.im = im;
		}
		
		// method to return the real part
		public double getReal() {return this.r; }
		public double getImaginary () {return this.im; }
		
		@Override
		public String toString () { return r + (im > 0.0 ? " + i" : " - i") + Math.abs(im); }
	}
	
	// EXCEPTION CLASSES
	private class PoincareException extends Exception
	{
		public PoincareException (String message) {
			super (message);
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		double x0, xf, y0, yf;
		x0 = -3.0;
		xf = 3.0;
		y0 = -3.0;
		yf = 3.0;
		
		// TODO Auto-generated method stub
		javax.swing.JFrame frame = new javax.swing.JFrame (
				"Poincare Crossing Program for BVP Oscillator");
		
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(new PoincareCrossing(frame, x0, xf, y0, yf));
		frame.setSize (580, 570);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
	}

}