package JosephsonJunction;

import java.awt.Color;
import java.awt.Graphics;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.JFrame;
import javax.swing.JPanel;

/* JJDynamics
 * This software is used to understudy the behaviour of 
 * a ac biased Josephson Junction and to investigate the system
 * at the chaotic state
 */

public class ResistiveJJ extends JPanel {
	
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
	private double alpha = beta;
	private static double omega = 0.60;
	private double i0 = 1.20;
	
	private static double pi = Math.PI;
	private static double T = pi / omega;
	private double dt = 0.005;
	
	private int n = 15000;
	private double [] x1, x2, y1, y2;
	
	// constructor
	public ResistiveJJ (double x0, double xf, double y0, double yf)
	{
		this.x0 = x0;
		this.xf = xf;
		this.y0 = y0;
		this.yf = yf;
		
		x1 = new double [n + 1];
		x2 = new double [n + 1];
		y1 = new double [n + 1];
		y2 = new double [n + 1];
		
		// start a new thread to solve iterates
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {
			
			public void run () {
				// set the initial values
				double v [] = new double [] {T / 5, 1.5, T / 10, 0.5};
				x1 [0] = v [0];
				y1 [0] = v [1];
				x2 [0] = v [2];
				y2 [0] = v [3];
				
				// solve for nth iterates
				for (int i = 0; i < n; ++i) {
					v = rungeKutta (v, i);
					
					// assign x and y
					x1 [i + 1] = v [0];
					y1 [i + 1] = v [1];
					x2 [i + 1] = v [2];
					y2 [i + 1] = v [3];
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
		PrintWriter writer = null;
		try {
			writer = new PrintWriter (new BufferedWriter (
					new FileWriter ("C:/Users/Olugbenga/Desktop/jj.txt")));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		writer.write ("y1\ty2\n");
		
		// plot the points
		for (int i = 1; i < x1.length; ++i) {
			int xx1, xx2, yy1, yy2;
			xx1 = toXPixels (i - 1);
			xx2 = toXPixels (i);
			yy1 = toYPixels (y1 [i - 1] - y2 [i - 1]);
			yy2 = toYPixels (y1 [i] - y2 [i]);
			
			g.drawLine (xx1, yy1, xx2, yy2);
			writer.write("\n" + y1 [i] + "\t" + y2 [i]);
			
		}
		writer.flush();
	}
	
	// method to return the x-coordinates
	private int toXPixels (int i)
	{
		return (xStart + (int) (i * width / x1.length) );
	}
	// method to convert the y-coordinates to pixels
	private int toYPixels (double y)
	{
		return (yStart + (int) (height * (y - y0) / (yf - y0)));
	}	// end method toPixels
	
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
				Math.sin(v [0]) - alpha * (v [1] - v [3]);
		vv [2] = v [3];
		vv [3] = i0 * Math.cos(omega * t) - beta * v [3] -
				Math.sin(v [2]) - alpha * (v [3] - v [2]);
		return vv;
	}
	
	// launch application
	public static void main(String[] args) {
		// instantiate the Frame to draw the figures
		ResistiveJJ app = new ResistiveJJ (0, 75, -300, 300);
		
		JFrame frame = new JFrame ("Dynamics of a Resistively" +
				" Coupled Josephson Junction");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setSize(500, 520);
		frame.add(app);
		frame.setVisible(true);
	}
}