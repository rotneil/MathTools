package phDModels;

import java.awt.Color;
import java.awt.Graphics;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.JFrame;
import javax.swing.JPanel;

/* JJDynamics
 * This software is used to understudy the behaviour of 
 * a ac biased Josephson Junction and to investigate the system
 * at the chaotic state
 */

public class JJDynamics extends JPanel {
	
	private static final long serialVersionUID = 1L;

	// instance variables
	private static int xStart = 30;
	private static int yStart = 15;
	private static int width = 430;
	private static int height = 430;
	
	// boundary parameters
	private double x0, xf, y0, yf;
	
	// parameter of equation
	private double beta = 0.25;
	private static double omega = 0.60;
	private double i0 = 1.20;
	
	private static double pi = Math.PI;
	private static double T = pi / omega;
	private double dt = 0.005;
	
	private int n = 150000;
	private double [] x, y;
	
	// constructor
	public JJDynamics (double x0, double xf, double y0, double yf)
	{
		this.x0 = x0;
		this.xf = xf;
		this.y0 = y0;
		this.yf = yf;
		
		x = new double [n + 1];
		y = new double [n + 1];
		
		// start a new thread to solve iterates
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {
			
			public void run () {
				// set the initial values
				double v [] = new double [] {T / 5, 1.5};
				x [0] = v [0];
				y [0] = v [1];
				
				// solve for nth iterates
				for (int i = 0; i < n; ++i) {
					v = rungeKutta (v, i);
					
					// assign x and y
					x [i + 1] = v [0];
					y [i + 1] = v [1];
					
					// adjust x
					if (Math.abs(x [i + 1]) > T) {
						int np = (int) (x [i + 1] / (2 * T) + Math.signum (x [i + 1]) * 0.5);
						x [i + 1] -= np * 2 * T;
					}
				}
				// show the series
				repaint ();
			}
		});
		
	}	// end of constructor
	
	// method to paint component
	public void paintComponent (Graphics g)
	{
		super.paintComponent(g);
		
		// draw the frame
		g.setColor(Color.white);
		g.fillRect(xStart, yStart, width, height);
		g.setColor(Color.black);
		g.drawRect(xStart, yStart, width, height);
		
		// show the legend
		for (int i = 1; i < 10; ++i) {
			g.drawString("_", xStart - 2, 
				yStart + (int)(i * height / 10));
			g.drawString("|", 
				xStart + (int)(i * width / 10), yStart + height + 3);
		}
		// show the values
		g.drawString("" + x0, (xStart - 3), yStart + height + 15);
		g.drawString("" + xf, xStart + width - 13, yStart + height + 15);
		g.drawString("" + y0, xStart - 25, yStart + height);
		g.drawString("" + yf, xStart - 25, yStart + 10);
		
		// plot the points
		for (int i = 0; i < x.length; ++i) {
			int px = (int) (xStart + width * (x [i] - x0) / (xf - x0));
			int py = (int) (yStart + height * (y [i] - y0) / (yf - y0));
			
			g.drawOval(px, py, 1, 1);
		}
	}
	
	// method rungeKutta
	private double [] rungeKutta (double [] v, int iterate)
	{
		// define the intermittent steps
		int l = v.length;
		double t = iterate * dt;
		double [] k1 = new double [l];
		double [] k2 = new double [l];
		double [] k3 = new double [l];
		double [] k4 = new double [l];
		
		// evaluate the steps
		k1 = f (v, t);
		
		for (int i = 0; i < l; ++i) 
			k2 [i] = v [i] + dt * k1 [i] / 2;
		k2 = f (k2, t + dt / 2);
		
		for (int i = 0; i < l; ++i)
			k3 [i] = v [i] + dt * k2 [i] / 2;
		k3 = f (k3, t + dt / 2);
		
		for (int i = 0; i < l; ++i)
			k4 [i] = v [i] + dt * k3 [i];
		k4 = f (k4, t + dt);
		
		// store the new value in k4
		for (int i = 0; i < l; ++i)
			k4 [i] = v [i] + dt * 
				(k1 [i] + 2 * (k2 [i] + k3 [i]) + k4 [i] ) / 6.0;
		return k4;
	}
	
	// method function
	private double [] f (double [] v, double t)
	{
		// define the new range
		double vv [] = new double [v.length];
		vv [0] = v [1];
		vv [1] = i0 * Math.cos(omega * t) - beta * v [1] -
				Math.sin(v [0]);
		return vv;
	}
	
	// launch application
	public static void main(String[] args) {
		// instantiate the Frame to draw the figures
		JJDynamics app = new JJDynamics (-T, T, -8, 8);
		
		JFrame frame = new JFrame ("Dynamics of a normalized" +
				" Josephson Junction");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setSize(500, 520);
		frame.add(app);
		frame.setVisible(true);
	}
}