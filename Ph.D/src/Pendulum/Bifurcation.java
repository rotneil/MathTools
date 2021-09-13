package Pendulum;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


public class Bifurcation extends javax.swing.JPanel 
	implements MouseMotionListener, MouseListener {
	
	// instance variables
	private int xStart = 50;
	private int yStart = 20;
	private int height = 450;
	private int width = 450;
	
	private static double x0 = 0.0;
	private static double xf = 2.78;
	private static double y0 = -Math.PI;
	private static double yf = Math.PI;
	
	private int grid = 1000;
	
	private double [] [] xy;
	private double rho, c = 0.05;
	private double dt = 0.01;
	
	// variables used by mouse event
	private int mx1, mx2, my1, my2;
	private boolean draw = false;
	
	public Bifurcation (javax.swing.JFrame frame, final double x0, final double xf, 
			final double y0, final double yf)
	{
		frame.setTitle("The bifurcation Diagram of Forced Damped Pedulum where c = " + c + 
				" for rho in [" + x0 + ", " + xf + "] and x in [" + y0 + ", " + yf + "]");
		this.x0 = x0;
		this.xf = xf;
		this.y0 = y0;
		this.yf = yf;
		
		xy = new double [grid] [grid];
		
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {

			public void run() 
			{
				// iterate through the horizontal values of rho
				for (int i = 0; i < grid; ++i) 
				{
					// assign value to rho
					rho = x0 + i * (xf - x0) / grid;
					
					// iterate through the vertical values of x
					for (int j = 0; j < grid; ++j) {
						// choose random initial condition
						double x = y0 + j * (yf - y0) / grid;
						double y = -5.0 + 10.0 * Math.random();
						
						// call function to give the steady state value of point
						xy [i] [j] = poincare (new double [] {x, y});
					}
					repaint ();
				}	// end for loop
				
				repaint ();
			}	// end method run
			
		});	// end executor
		
		// add mouse motion listener
		addMouseListener (this);
		addMouseMotionListener (this);
	}
	
	// method to state the steady state value of trajectory
	private double poincare (double [] v)
	{
		int iterate = 50000;
		// perform n-iterate of point v
		for (int i = 0; i < iterate; ++i){
			v = rungeKutta (v, dt * (i + 1));
			
			// always bring x back
			v [0] = torus (v[0], 2.0 * Math.PI);
		}
		
		// iterate through to get the first crossing of the Poincare section
		double [] v2, vv;
		for (int i = iterate; i < iterate + 2000; ++i) {
			v = rungeKutta (v, dt * i);
			v [0] = torus (v[0], 2.0 * Math.PI);
			
			// compute the next iterate
			v2 = rungeKutta (v, dt * (i + 1));
			v2[0] = torus (v2[0], 2.0 * Math.PI);
			
			vv = f (v, dt * i);
			if (vv[1] < 0 && v[1] > 0 && v2[1] <= 0)
				return (v[0] + v2[0]) / 2.0;
		}
		return 2.0 * Math.PI;
	}
	
	// equation definition
	private double[] f (double [] v, double t)
	{
		double [] vv = new double [v.length];
		
		vv[0] = v [1];
		vv[1] = -c * v[1] - Math.sin(v[0]) + rho * Math.sin(t);
		
		return vv;
	}
	
	// bring x and y-values back to the region [-pi, pi]
	private double torus (double z, double mod)
	{
		int sign = z >= 0 ? 1 : -1;
		int np = (int) (z / mod + sign * 0.5);
		z -= mod * np;
		
		return z;
	}
	
	// Runge Kutta fourth order numerical integeration
	private double [] rungeKutta (double [] v, double t)
	{
		int l = v.length;
		double [] c1, c2, c3, c4;
		
		// initialize the intermediate steps
		c1 = new double [l];
		c2 = new double [l];
		c3 = new double [l];
		c4 = new double [l];
		
		c1 = f(v, t);
		
		for (int i = 0; i < l; ++i)
			c2[i] = v[i] + dt * c1[i] / 2;
		c2 = f(c2, t + dt / 2);
		
		for (int i = 0; i < l; ++i)
			c3[i] = v[i] + dt * c2[i] / 2;
		c3 = f(c3, t + dt / 2);
		
		for (int i = 0; i < l; ++i)
			c4[i] = v[i] + dt * c3[i];
		c4 = f(c4, t + dt);
		
		for (int i = 0; i < l; ++i)
			c1[i] = v[i] + dt * (c1[i] + 2 * (c2[i] + c3[i]) + c4[i]) / 6.0;
		
		return c1;
	}
	
	// method to paint on the panel
	public void paintComponent (java.awt.Graphics g)
	{
		super.paintComponent(g);
		
		// make the draw area white
		setBackground (java.awt.Color.white);
		g.drawRect(xStart, yStart, width, height);
		
		// show the legend
		g.setColor(java.awt.Color.black);
		g.drawString("pi", xStart - 25, yStart + 5);
		g.drawString("-pi", xStart - 25, yStart + height + 5);
		g.drawString("" + x0, xStart - 5, yStart + height + 15);
		g.drawString("" + xf, xStart + width - 10, yStart + height + 15);
		
		// show minor marks
		for (int i = 0; i < 21; ++i) {
			// divide y-axis
			g.drawString("-", xStart, yStart + (int) (i * height / 20) + 5);
			g.drawString("-", xStart + width - 2, yStart + (int) (i * height / 20) + 5);
			g.drawString("'", xStart + (int) (i * width / 20), yStart + 10);
			g.drawString("'", xStart + (int) (i * width / 20), yStart + height + 7);
		}
		
		// draw plot on the panel
		for (int i = 0; i < xy.length; ++i) {
			int xx = xStart + (int) (i * width / xy.length);
			for (int j = 0; j < xy [i].length; ++j) {
				int yy = yStart + (int) (
						height * (yf - xy [i] [j]) / (yf - y0) );
				
				// draw a point here
				g.drawOval(xx, yy, 1, 1);
			}
		}
		
		if (draw)
			g.drawRect(Math.min(mx1, mx2), Math.min(my1, my2),
					Math.abs(mx1 - mx2), Math.abs(my1 - my2));
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		javax.swing.JFrame frame = new javax.swing.JFrame ();
		//frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(new Bifurcation (frame, x0, xf, y0, yf));
		frame.setSize (580, 560);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
	}

	public void mouseClicked(MouseEvent e) {
		// TODO Auto-generated method stub
		draw = false;
	}

	public void mousePressed(MouseEvent e) {
		// TODO Auto-generated method stub
		mx1 = e.getX();
		my1 = e.getY();
		draw = false;
	}

	// method to give a 2 s.f double
	private double to4sf (double x)
	{
		double sf = Math.floor(10000 * x + 0.5);
		return sf / 10000.0;
	}
	
	public void mouseReleased(MouseEvent e) {
		// TODO Auto-generated method stub
		// start another frame
		if (draw)
			startAnotherFrame (
				getXValue (Math.min(mx1, mx2)),
				getXValue (Math.max(mx1, mx2)),
				getYValue (Math.max(my1, my2)),
				getYValue (Math.min(my1, my2)));
	}

	public void mouseEntered(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	public void mouseExited(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	public void mouseDragged(MouseEvent e) {
		// TODO Auto-generated method stub
		mx2 = e.getX();
		my2 = e.getY();
		draw = true;
	}

	public void mouseMoved(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	private double getXValue (int xx)
	{
		return (x0 + (xf - x0) * (xx - xStart) / width);
	}
	
	private double getYValue (int yy)
	{
		return (yf + (y0 - yf) * (yy - yStart) / height);
	}
	
	private void startAnotherFrame (
			double x0, double xf, double y0, double yf)
	{
		x0 = to4sf(x0);
		xf = to4sf(xf);
		y0 = to4sf(y0);
		yf = to4sf(yf);
		
		// TODO Auto-generated method stub
		javax.swing.JFrame frame = new javax.swing.JFrame ();
		frame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
		frame.add(new Bifurcation (frame, x0, xf, y0, yf));
		frame.setSize (580, 560);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
	}
}