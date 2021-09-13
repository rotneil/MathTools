package bifurcation;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.JFrame;

public class StabilityTestII extends javax.swing.JPanel
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
	
	// constructor
	public StabilityTestII (javax.swing.JFrame frame,
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
		
		// add mouse click listener
		addMouseListener(new java.awt.event.MouseAdapter () {
			public void mouseClicked (java.awt.event.MouseEvent evt) {
				final double k = toY (evt.getY());
				final double g = toX (evt.getX());
				
				ExecutorService executor = Executors.newCachedThreadPool();
				executor.execute(new Runnable () {
					public void run () {
						solve (g, k);
					}
				});
				executor.shutdown();
			}
		});
	}
	
	// method to perform computation
	private void solve (double g, double k)
	{
		System.out.printf("\n[g=%f,k=%f]", g, k);
		
		// set the initial condition
		double [] v = new double [] {
				3.0 * Math.random(), 3.0 * Math.random(), 1.0, 0.0, 0.0, 1.0};
		
		// check for non-periodic orbit
		if (isNonPeriodic (g, k, v)) {
			System.out.println ("\nNon periodic");
			return;
		}
		
		// get the fixed on Poincare section
		Point fp = getPoincareFixedPoint (g, k, v);
		
		// reset the variational coeficients and integrate until the next crossing
		v[0] += 0.01; v[2] = 1.0; v[3] = 0.0; v[4] = 0.0; v[5] = 1.0;
		v = getNextPoincareCrossing (g, k, v);
		
		// get the eigen values of the variational matrix
		String [] h = getEigenValues (new double [] [] {{v[2], v[3]}, {v[4], v[5]}});
		
		System.out.printf("\n%f\th1 = %s\th2 = %s", v[0], h[0], h[1]);
	}
	
	// method to check if an orbit is periodic
	private boolean isNonPeriodic (double g, double k, double [] v)
	{
		// get the fixed point
		java.util.ArrayList<Point> fixedPt = getFixedPoint(g, k);
		
		// iterate for 10000 times
		for (int j = 0; j < 80000; ++j)
			v = rungeKutta (g, k, v);
		
		// confirm that orbit is not attracted to the fixed points
		for (Point pt : fixedPt)
			if (Math.abs(pt.getX() - v[0]) <= 0.0001 && Math.abs(pt.getY() - v[1]) <= 0.0001)
				return true;
		
		return false;
	}
	
	// method to find and return the fixed point on the Poincare section
	private Point getPoincareFixedPoint (double g, double k, double [] v)
	{
		double [] vv = new double [v.length];
		
		while (true) {
			vv = getNextPoincareCrossing (g, k, v);
			
			if (Math.abs(v[0] - vv[0]) <= 0.0001)
				return new Point (vv[0], vv[1]);
			
			v = vv;
		}
	}
	
	/**
	 * This method searches for the next Poincare section crossing
	 * @param theta
	 * @param v
	 * @return
	 */
	private double [] getNextPoincareCrossing (double g, double k, double [] v)
	{	
		// iterate until there is a crossing
		while (true) {
			v = rungeKutta (g, k, v);
			
			// get the direction vector
			double [] vt = f(g, k, v);
			if (vt[1] > 0.0 && Math.abs(v[1]) <= 0.0001)
				return v;
		}
	}
    
	// method to set the Poincare sectio angle
	private double getPoincareAngle (double g, double k)
	{
		// get the fixed points
		java.util.ArrayList<Point> fixedPoints = getFixedPoint (g, k);
		if (fixedPoints.size() == 1)
			return 0.0;
		else {
			double x1, x2, y1, y2;
			x1 = fixedPoints.get(1).getX();
			x2 = fixedPoints.get(2).getX();
			y1 = fixedPoints.get(1).getY();
			y2 = fixedPoints.get(2).getY();
			
			return Math.atan((y2 - y1) / (x2 - x1));
		}
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
	 		fixedPoint.add(new Point (m, n));
	 		fixedPoint.add(new Point (-m, -n));
        }
        
 		return fixedPoint;
    }
	
    // Recursive Newton-Raphson method call to 
 	private double refinePoint (double g, double k, double x0)
 	{
 		// local variable
 		double fx = Math.tanh(g * x0) - x0 / k;
 		double ffx = g * Math.pow(Math.cosh(g * x0), -2.0) - 1.0 / k;
 		double x1 = x0 - fx / ffx;
 		
 		if (Math.abs(x1 - x0) <= 0.0001)
 			return x1;
 		else
 			return refinePoint (g, k, x1);
 	}
 	
 	// method to get the eigen values from 
	private String [] getEigenValues (double [] [] a)
	{
		// local variable
		java.text.DecimalFormat df = new java.text.DecimalFormat("0.0000;0.0000");
		double r1, r2, i;
		double tr, det;
		tr = a[0][0] + a[1][1];
		det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
		double d = sq (tr) - 4.0 * det;
		
		if (d >= 0.0) {
			r1 = 0.5 * (tr + Math.sqrt(d));
			r2 = 0.5 * (tr - Math.sqrt(d));
			return new String [] {"" + df.format (r1), "" + df.format (r2)};
		} else {
			r1 = 0.5 * tr;
			i = 0.5 * Math.sqrt(-d);
			return new String [] {df.format (r1) + " + i" + df.format (i), 
					df.format (r1) + " - i" + df.format (i)};
		}
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
		return (v[0] * vt[1] - v[1] * vt[0]) / (v[0] * v[0] + v[1] * v[1]);
	}
	
	// method to return the radius
	private double getRadius (double [] v) 
	{
		return Math.sqrt(v[0] * v[0] + v[1] * v[1]);
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
	
	private double sq (double x) {return x * x; }
	
	// method to show result in certain decimal places
	private double todp (double x, int p)
	{
		x *= Math.pow(10, p) + 0.5;
		int value = (int) x;
		return (1.0 * value / Math.pow(10, p));
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
		javax.swing.JFrame frame = new javax.swing.JFrame ("Stability Test Program");
		
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(new StabilityTestII (frame, x0, xf, y0, yf));
		frame.setSize (580, 560);;
		frame.setVisible(true);
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
}
