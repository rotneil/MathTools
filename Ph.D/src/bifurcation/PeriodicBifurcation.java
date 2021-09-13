package bifurcation;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.JFrame;

public class PeriodicBifurcation extends javax.swing.JPanel 
{
	// instance variables
	private javax.swing.JFrame mFrame;
	
	private int xStart = 55;
	private int yStart = 20;
	private int height = 450;
	private int width = 450;
	
	// plot dimension
	private final double x0;
	private final double xf;
	private final double y0;
	private final double yf;
	private final int grid = 1000;
	
	private final short [] [] xy;
	private final double dt = 0.01;
	private final double epsilon = 0.0001;
	private final double deltaX = 0.01;
	
	private final short NO_BIFURCATION = 0;
	private final short SADDLE_NODE_BIFURCATION = 1;
	private final short PERIOD_DOUBLING_BIFURCATION = 2;
	private final short PITCHFORK_BIFURCATION = 3;
	private final short TRANSCRITICAL_BIFURCATION = 4;
	private final short HOPF_BIFURCATION = 5;
	private final short NEIMARK_SACKER_BIFURCATION = 6;
	
	
	private String status = "Status: Started ...";
	
	// constructor
	public PeriodicBifurcation (javax.swing.JFrame frame,
			final double x0, final double xf, final double y0, final double yf)
	{
		// instantiate the instance data
		this.mFrame = frame;
		frame.setTitle("Two Parameter Bifurcation diagram for BVP Oscillator " +
				"with dx = " + deltaX + " and e = " + epsilon);
		setBackground (java.awt.Color.white);
		
		this.x0 = x0;
		this.xf = xf;
		this.y0 = y0;
		this.yf = yf;
		
		xy = new short [grid + 1] [grid + 1];
		
		// start computing the bifurcation patterns
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {
			public void run () {
				// space variable
				double g, k;
				
				// loop through the whole space
				for (int i = 500; i < xy.length; ++i) {
					// assign the horizontal value
					g = x0 + i * (xf - x0) / grid;
					
					// iterate along the vertical axis
					for (int j = 250; j <= 500; ++j) {//xy[i].length; ++j) {
						k = y0 + j * (yf - y0) / grid;
						
						// compute the bifurcation type for this instance
						try {
							xy[i][j] = getBifurcationType (g, k);
						} catch (PoincareException e) {}
					
						// compute the status of computation and Level
						status = "Status: " + (int) (100.0 * i / grid) + "%   " +
								"Step: " + (i % 10) + ":   " + 
								(int) (100.0 * j / grid) + "%";//grid) + "%";
						
						// show the basin
						updatePaint ();
					}
					
					// compute the status of computation and Level
					status = "Status: " + (int) (100.0 * i / grid) + "%   " +
							"Step: " + (i % 10) + ":   100%";
					
					// show the basin
					updatePaint ();
				}	// end outter for-loop
				updatePaint ();
			}
		});
		executor.shutdown();
	}
	
	/**
	 * This method detects the bifurcation thresholds for all the 
	 * bifurcation types
	 * @param g
	 * @param k
	 * @return
	 */
	private short getBifurcationType (double g, double k) throws PoincareException
	{
		// assign initial value for the orbit
		double [] v = new double [] {5.5, 3.8, 1.0, 0.0, 0.0, 1.0};
		
		// get Poincare section
		double section = getPoincareSection (g, k);
		
		// get Poincare fixed point
		double [] vp = getPoincareFixedPoint (section, g, k, v);
		
		// get return time
		int returnTime = getReturnTime (section, g, k, vp);
		
		// set new orbit with Ne(e)
		// vp[0] += 0.1; vp[1] += 0.1 * tan (section);
		vp[2] = 1.0; vp[3] = 0.0; vp[4] = 0.0; vp[5] = 1.0;
		
		// iterate the new orbit for returnTime-times
		for (int counter = 0; counter <= returnTime; ++ counter)
			vp = rungeKutta (g, k, vp);
		
		// compute the eigen values and test for bifurcatio thresholds
		for (ComplexNumber h: getEigenValues (vp)) {
			double l = h.getModulus();
			double Re = h.getReal();
			double Im = h.getImaginary();
			
			// Andronov Hopf bifurcation
			if (Math.abs(Re) <= 0.01 && Im != 0.0)
				return HOPF_BIFURCATION;
			if (Math.abs(Re + 1.0) <= 0.01 && Math.abs(Im) <= 0.01)
				return PERIOD_DOUBLING_BIFURCATION;
			if (Math.abs(Re - 1.0) <= 0.01 && Math.abs(Im) <= 0.01)
				return SADDLE_NODE_BIFURCATION;
			if (Math.abs(l - 1.0) <= 0.01)
				return NEIMARK_SACKER_BIFURCATION;
		}
		
		return NO_BIFURCATION;
	}
	
	/**
	 * This method iterates an orbit through one cycle of the limit curve
	 * and returns the time (in iterations) required to make such a cycle
	 * @param section
	 * @param g
	 * @param k
	 * @param v
	 * @return
	 */
	private int getReturnTime (double section, double g, double k, double [] v0)
	{
		// initial values
		double [] v = new double [v0.length];
		double [] vv = new double [v0.length];
		v = rungeKutta (g, k, v0);
		int period = 1;
		
		// continuously iterate vv
		while (true) {
			vv = rungeKutta (g, k, v);
			++period;
			
			// check for section crossing
			if (v[1] < v[0] * Math.tan(section) && vv[1] >= vv[0] * Math.tan(section))
				return period;
			
			v = vv;
		}
	}
	
	/**
	 * This method seeks for the fixed points on Poincare Section
	 * @param section
	 * @param g
	 * @param k
	 * @param v
	 * @return
	 */
	private double [] getPoincareFixedPoint (double section, double g, double k, double [] v0)
		throws PoincareException
	{
		// local variable
		double [] v1 = new double [v0.length];
		double [] v2 = new double [v0.length];
		double [] v3 = new double [v0.length];
		double [] v4 = new double [v0.length];
		double [] v5 = new double [v0.length];
		double [] v6 = new double [v0.length];
		
		v1 = getNextPoincareCrossing (section, g, k, v0);
		v2 = getNextPoincareCrossing (section, g, k, v1);
		v3 = getNextPoincareCrossing (section, g, k, v2);
		v4 = getNextPoincareCrossing (section, g, k, v2);
		v5 = getNextPoincareCrossing (section, g, k, v2);
		
		// continously iterate orbit until a fixed point is found on the Poincare
		while (true) {
			v6 = getNextPoincareCrossing (section, g, k, v5);
			
			// check for fixed point
			if (Math.abs(v2[0] - v1[0]) <= 0.001 && Math.abs(v2[1] - v1[1]) <= 0.001)
				return v2;
			if (Math.abs(v3[0] - v1[0]) <= 0.001 && Math.abs(v3[1] - v1[1]) <= 0.001)
				return v3;
			if (Math.abs(v4[0] - v1[0]) <= 0.001 && Math.abs(v4[1] - v1[1]) <= 0.001)
				return v4;
			if (Math.abs(v5[0] - v1[0]) <= 0.001 && Math.abs(v5[1] - v1[1]) <= 0.001)
				return v5;
			if (Math.abs(v6[0] - v1[0]) <= 0.001 && Math.abs(v6[1] - v1[1]) <= 0.001)
				return v6;
			
			// interchange orbit
			v1 = v2;
			v2 = v3;
			v3 = v4;
			v4 = v5;
			v5 = v6;
		}
	}
	
	/**
	 * This method iterates a orbit until it crosses the Poincare Section
	 * @param section
	 * @param g
	 * @param k
	 * @param v
	 * @return
	 */
	private double [] getNextPoincareCrossing (double section, double g, double k, double [] v0)
		throws PoincareException
	{
		// initial values
		double [] v = new double [v0.length];
		double [] vv = new double [v0.length];
		java.util.ArrayList<Point> fixedPoint = getFixedPoint(g, k);
		
		v = rungeKutta (g, k, v0);
		
		// continuously iterate vv
		while (true) {
			vv = rungeKutta (g, k, v);
			
			// check for section crossing
			if (v[1] < v[0] * Math.tan(section) && vv[1] >= vv[0] * Math.tan(section))
				return vv;
			else
				for (Point p : fixedPoint)
					if (Math.abs(vv[0] - p.getX()) <= 0.0000001 && 
							Math.abs(vv[1] - p.getY()) <= 0.0000001)
						throw new PoincareException ("No Poincare Section crossing");
			
			v = vv;
		}
	}
	
	/**
	 * This method returns the angle of the Poincare section used
	 * for this value of g and k
	 * @param g
	 * @param k
	 * @return
	 */
	private double getPoincareSection (double g, double k)
	{
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
    
    /**
     * This method computes the eigen values of a matrix
     * @param v
     * @return
     */
    private ComplexNumber [] getEigenValues (double [] v)
    {
    	// local variables
    	double m1 = v[2], m2 = v[3], m3 = v[4], m4 = v[5];
    	double tr = m1 + m4;
    	double det = m1 * m4 - m2 * m3;
    	double d = sq(tr) - 4 * det;
    	
    	if (d >= 0)
    		return new ComplexNumber [] {
				new ComplexNumber (0.5 * (tr + Math.sqrt(d)), 0.0),
				new ComplexNumber (0.5 * (tr - Math.sqrt(d)), 0.0) };
    	else
    		return new ComplexNumber [] {
    			new ComplexNumber (0.5 * tr, 0.5 * Math.sqrt(-d)),
    			new ComplexNumber (0.5 * tr, -0.5 * Math.sqrt(-d))};
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
 		
 		// variational components
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
				case NEIMARK_SACKER_BIFURCATION:
					g.setColor(java.awt.Color.orange);
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
		public double getModulus () {return Math.sqrt(sq(r) + sq (im)); }
		
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
		frame.add(new PeriodicBifurcation (frame, x0, xf, y0, yf));
		frame.setSize (580, 560);
		//frame.setLocationRelativeTo(null);
		frame.setVisible(true);
	}

}