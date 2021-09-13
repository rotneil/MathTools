package bifurcation;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

public class TwoParameterBifurcation extends javax.swing.JPanel
{
	// instance variables
	private javax.swing.JFrame mFrame;
	
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
	
	private short [] [] xy;
	private double dt = 0.01;
	private double e = 0.001;
	
	private final short NO_BIFURCATION = 0;
	private final short SADDLE_NODE_BIFURCATION = 1;
	private final short PERIOD_DOUBLING_BIFURCATION = 2;
	private final short PITCHFORK_BIFURCATION = 3;
	private final short TRANSCRITICAL_BIFURCATION = 4;
	private final short HOPF_BIFURCATION = 5;
	
	private String status = "";
	
	// constructor
	public TwoParameterBifurcation (javax.swing.JFrame frame,
			final double x0, final double xf, final double y0, final double yf)
	{
		// instantiate the instance data
		this.mFrame = frame;
		setBackground (java.awt.Color.white);
		
		this.x0 = x0;
		this.xf = xf;
		this.y0 = y0;
		this.yf = yf;
		
		xy = new short [grid] [grid];
		
		// start computing the bifurcation patterns
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {
			public void run () {
				// declare the parameters
				double g, k;
				
				// perform iteration
				for (int i = 500; i < grid; ++i) 
				{
					// assign the horizontal value
					g = x0 + i * (xf - x0) / grid;
					
					// assign the vertical value
					for (int j = 250; j < 500; ++j) 
					{
						k = y0 + j * (yf - y0) / grid;
						
						// set an initial condition
						double [] v = new double [] {2.5, 2.4, 1.0, 0.0, 0.0, 1.0};
						
						// set the region
						try {
							if (g * k > 1.0)
								xy[i][j] = getBifurcationPoint (g, k, v);
							else
								xy [i] [j] = NO_BIFURCATION;
						} catch (PoincareException e) {
							System.out.println(e.getMessage() + 
									" at [g, k] = [" + g + ", " + k + "]");
						}
						
						// compute the status of computation and Level
						status = "Status: " + (int) (100.0 * (i - 500) / 500) + "%   " +
								"Step: " + (i % 5) + ":   " + 
								(int) (100.0 * (j - 250) / 250) + "%";
						
						// show the basin
						updatePaint ();
					}
					
					// compute the status of computation and Level
					status = "Status: " + (int) (100.0 * (i - 500) / 500) + "%   " +
							"Step: " + (i % 5) + ":   100%";
					
					// show the basin
					updatePaint ();
				}	// end for loop
				
				updatePaint ();
			}
		});
		executor.shutdown();
		
	}
	
	private void updatePaint () 
	{
		SwingUtilities.invokeLater(new Runnable () {
			public void run () {
				repaint();
			}
		});
	}
	
	// method to perform computation
	private short getBifurcationPoint (double g, double k, double [] v) throws PoincareException
	{
		// get the poincare section angle
		double theta = getPoincareSectionAngle (g, k);
		
		// get the fixed on Poincare section and the return time
		v = getPoincareFixedPoint (theta, g, k, v);
		int period = getReturnTime (theta, g, k, v);
		
		// reset the variational coeficients and integrate until the next crossing
		v[0] += 0.01; v[1] += 0.01 * Math.tan(theta);
		v[2] = 1.0; v[3] = 0.0; v[4] = 0.0; v[5] = 1.0;
		
		// iterate the orbit until it crosses Poincare section
		for (int i = 0; i < period; ++i)
			v = rungeKutta (g, k, v);
		
		// get the eigen values of the variational matrix
		ComplexNumber [] eigenValues = getEigenValues (new double [] [] {{v[2], v[3]}, {v[4], v[5]}});
		
		// check for bifurcation threshold
		for (ComplexNumber h: eigenValues)
			if (Math.abs(h.getReal() - 1) <= 0.0001 && Math.abs(h.getImaginary()) <= 0.0001)
				return SADDLE_NODE_BIFURCATION;
			else if (Math.abs(h.getReal() + 1) <= 0.0001 && Math.abs(h.getImaginary()) <= 0.0001)
				return PERIOD_DOUBLING_BIFURCATION;
			else if (Math.abs(h.getReal()) <= 0.0001 && h.getImaginary() != 0.0)
				return HOPF_BIFURCATION;
		
		return NO_BIFURCATION;
	}
	
	// this method returns the return time
	private int getReturnTime (double theta, double g, double k, double [] v)
		throws PoincareException
	{
		int period = 0;
		double [] vv = new double [v.length];
		for (int i = 0; i < v.length; ++i)
			vv[i] = v[i];
		
		while (true) {
			vv = rungeKutta (g, k, vv);
			++period;
			
			// get the direction vector and then check for Poincare crossing
			double [] vt = f(g, k, vv);
			if (getAngleSpeed (vv, vt) > 0.0 && Math.abs(theta - getAngle (vv)) <= 0.0001)
				return period;
			// search for fixed point basin
			else
				for (Point p : getFixedPoint (g, k))
					if (Math.abs(p.getX() - vv[0]) <= 0.000001 && 
						Math.abs(p.getY() - vv[1]) <= 0.000001)
						throw new PoincareException ("A Fixed point detected at " + p.toString());
		}
	}
	
	// method to find and return the fixed point on the Poincare section
	private double [] getPoincareFixedPoint (double theta, double g, double k, double [] v)
		throws PoincareException
	{
		double [] v0 = new double [v.length];
		double [] vv = new double [v.length];
		
		// assign the initial crossing
		v0 = getNextPoincareCrossing (theta, g, k, v);
		
		// iterate crossings
		while (true) {
			vv = getNextPoincareCrossing (theta, g, k, v0);
			
			if (Math.abs(v0[0] - vv[0]) <= 0.0001 && Math.abs(v0[1] - vv[1]) <= 0.0001)
				return vv;
			
			v0 = vv;
		}
	}
	
	/**
	 * This method searches for the next Poincare section crossing
	 * @param theta
	 * @param v
	 * @return
	 * @throws PoincareException if orbit falls to a fixed point
	 */
	private double [] getNextPoincareCrossing (double theta, double g, double k, double [] v)
		throws PoincareException
	{
		double [] vv = new double [v.length];
		for (int i = 0; i < v.length; ++i)
			vv[i] = v[i];
		
		// iterate until there is a crossing
		while (true) {
			vv = rungeKutta (g, k, vv);
			
			// get the direction vector
			double [] vt = f(g, k, vv);
			if (getAngleSpeed (vv, vt) > 0.0 && Math.abs(theta - getAngle (vv)) <= 0.0001)
				return vv;
			// search for fixed point basin
			else
				for (Point p : getFixedPoint (g, k))
					if (Math.abs(p.getX() - vv[0]) <= 0.000001 && 
						Math.abs(p.getY() - vv[1]) <= 0.000001)
						throw new PoincareException ("A Fixed point detected at " + p.toString());
		}
	}
	
	// method to compute the best poincare section angle
	private double getPoincareSectionAngle (double g, double k)
	{
		java.util.ArrayList<Point> fixedPoint = getFixedPoint (g, k);
		
		if (fixedPoint.size() == 1)
			// return angle-45 in radians
			return 0.0;//Math.PI / 4.0;
	
		Point p1, p2;
		double x1, x2, y1, y2;
		p1 = fixedPoint.get(1);
		p2 = fixedPoint.get(2);
		x1 = p1.getX(); y1 = p1.getY();
		x2 = p2.getX(); y2 = p2.getY();
		
		return Math.atan((y2 - y1) / (x2 - x1));
	}
	
	// method to return the orbit angle in radian
	private double getAngle (double [] v)
	{
		if (v[0] == 0.0) {
			if (v[1] == 0.0)
				return 0.0;
			else
				return (v[1] > 0.0 ? Math.PI / 2.0 : 1.5 * Math.PI);
		} else if (v[0] > 0.0) {
			// first quadrant
			if (v[1] >= 0.0)
				return Math.atan(v[1] / v[0]);
			// fourth quadrant
			else
				return 2.0 * Math.PI - Math.atan(Math.abs(v[1] / v[0]));
		} else {
			// second quadrant
			if (v[1] >= 0.0)
				return Math.PI - Math.atan(Math.abs(v[1] / v[0]));
			// third quadrant
			else
				return Math.PI + Math.atan(Math.abs(v[1] / v[0]));
		}
	}
	
	// method to return the angle speed
	private double getAngleSpeed (double [] v, double [] vt)
	{
		return (v[0] * vt[1] - v[1] * vt[0]) / (sq(v[0]) + sq(v[1]));
	}
	
	// method to return the radius
	private double getRadius (double [] v) 
	{
		return Math.sqrt(sq(v[0]) + sq(v[1]));
	}
	
	// method to set fixed point
	private java.util.ArrayList<Point> getFixedPoint (double g, double k)
    {
    	// local variable
    	java.util.ArrayList<Point> fixedPoint = new java.util.ArrayList<Point>();
    	
    	// the origin is a constant fixed point
        fixedPoint.add(new Point ());
        double m = 0.0;
        
    	// the other points
        if (g * k > 1) {
        	m = Math.sqrt(3.0 * (g * k - 1.0) / (k * g * g * g));
	        m = refinePoint (g, k, m);
	 		
	        double n = m / k;
	 		fixedPoint.add(new Point (-m, -n));
	 		fixedPoint.add(new Point (m, n));
        }
        
 		return fixedPoint;
    }
	
	private double sq (double x) {return x * x; }
	
    // Recursive Newton-Raphson method call to 
 	private double refinePoint (double g, double k, double x0)
 	{
 		// local variable
 		if (g * k <= 1)
 			return 0.0;
 		
 		double fx = Math.tanh(g * x0) - x0 / k;
 		double ffx = g * Math.pow(Math.cosh(g * x0), -2.0) - 1.0 / k;
 		double x1 = x0 - fx / ffx;
 		
 		if (Math.abs(x1 - x0) <= 0.0001)
 			return x1;
 		else
 			return refinePoint (g, k, x1);
 	}
 	
 	// method to get the eigen values from 
	private ComplexNumber [] getEigenValues (double [] [] a)
	{
		// local variable
		double r1, r2, i;
		double tr, det;
		tr = a[0][0] + a[1][1];
		det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
		double d = sq (tr) - 4.0 * det;
		
		if (d >= 0.0) {
			r1 = 0.5 * (tr + Math.sqrt(d));
			r2 = 0.5 * (tr - Math.sqrt(d));
			return new ComplexNumber [] {new ComplexNumber (r1, 0.0), new ComplexNumber (r2, 0.0)};
		} else {
			r1 = 0.5 * tr;
			i = 0.5 * Math.sqrt(-d);
			return new ComplexNumber [] {new ComplexNumber (r1, i), new ComplexNumber (r1, -i)};
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
	
	private double [] f (double g, double k, double [] v)
	{
		// local variables
		double [] vv = new double [v.length];
		
		vv[0] = -v[1] + Math.tanh(g * v[0]);
		vv[1] = v[0] - k * v[1];
		
		vv[2] = g / sq(Math.cosh(g * v[0])) * v[2] - v[4];
		vv[4] = v[2] - k * v[4];
		
		vv[3] = g / sq(Math.cosh(g * v[0])) * v[3] - v[5];
		vv[5] = v[3] - k * v[5];
		
		return vv;
	}
	
	// method to paint on the panel
	public void paintComponent (java.awt.Graphics g)
	{
		super.paintComponent(g);
		
		// draw the bifurcation patterns on the screen
		g.drawLine(xStart, yStart + height, 
				(int)(xStart + width / 2.0), (int) (yStart + height / 2.0));
		int mesh = 100;
		for (int i = 1; i <= mesh; ++i) {
			double x1 = 0.5 + (i - 1) * 1.5 / mesh;
			double x2 = 0.5 + i * 1.5 / mesh;
			
			g.drawLine(toPx(x1), toPy(1.0 / x1), toPx(x2), toPy(1.0 / x2));
		}
			
		// make the draw area white
		g.drawRect(xStart, yStart, width, height);
		
		// show the bifurcation pattern
		for (int i = 0; i < xy.length; ++i) {
			int xx = xStart + (int) (i * width / grid);
			for (int j = 0; j < xy [i].length; ++j) {
				int yy = yStart + height - (int) (j * height / grid);
				
				// show the bifurcation pattern here
				switch (xy[i][j]) {
				case HOPF_BIFURCATION:
					g.setColor(java.awt.Color.blue);
					g.drawOval(xx, yy, 1, 1);
					break;
				case SADDLE_NODE_BIFURCATION:
					g.setColor(java.awt.Color.red);
					g.drawOval(xx, yy, 1, 1);
					break;
				case PITCHFORK_BIFURCATION:
					g.setColor(java.awt.Color.green);
					g.drawOval(xx, yy, 1, 1);
					break;
				default:
					g.setColor(java.awt.Color.white);
				}
			}
		}
		
		// show status
		g.setColor(java.awt.Color.black);
		g.drawString(status, 5, 515);
		
		// SHOW THE PLOT AREA
		g.drawRect(xStart, yStart, width, height);
		
		// show the legend
		g.drawString("" + yf, xStart - 25, yStart + 5);
		g.drawString("" + y0, xStart - 25, yStart + height + 5);
		g.drawString("" + x0, xStart - 5, yStart + height + 15);
		g.drawString("" + xf, xStart + width - 10, yStart + height + 15);
		
		// show minor marks
		for (int i = 0; i < 21; ++i) {
			// divide y-axis
			g.drawString("-", xStart, yStart + (int) (i * (height) / 20) + 4);
			g.drawString("-", xStart + width - 2, yStart + (int) (i * (height) / 20) + 4);
			g.drawString("'", xStart + (int) (i * width / 20) - 1, yStart + 10);
			g.drawString("'", xStart + (int) (i * width / 20) - 1, yStart + (height) + 7);
		}
	}
	
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
			return "[x=" + df.format(this.x) + ",y=" + df.format(this.y) + "]";
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
		x0 = 0.0;
		xf = 2.0;
		y0 = 0.0;
		yf = 2.0;
		
		// TODO Auto-generated method stub
		javax.swing.JFrame frame = new javax.swing.JFrame ("Two Parameter Bifurcation " +
				"diagram for BVP Oscillator");
		
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(new TwoParameterBifurcation (frame, x0, xf, y0, yf));
		frame.setSize (580, 560);
		//frame.setLocationRelativeTo(null);
		frame.setVisible(true);
	}

}
