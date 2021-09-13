package rossler;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Graphics;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;


public class RosslerBifurcation extends JPanel 
{
	// instance variables
	private static double y0 = 0.0;
	private static double yf = 70.0;
	private static double x0 = 1.0;
	private static double xf = 46.0;
	
	private int xStart = 50;
	private int yStart = 20;
	private int height = 500;
	private int width = 500;
	private int grid = 1000;
	private int iterate = 500;
	private static int screenWidth = 630;
	private static int screenHeight = 610;
	
	private double dt = 0.01;
	private double xy [][];
	
	// lorenz parameter
	private double c = x0;
	private double a = 0.1;
	private double b = 0.1;
	
	// Constructor
	public RosslerBifurcation ()
	{
		// instantiate the bif values
		xy = new double [grid + 1] [iterate];
		
		SwingUtilities.invokeLater(new Runnable () {
			public void run () {
				// iterate through the rho-values
				for (int i = 0; i < xy.length; ++i) {
					c = x0 + i * (xf - x0) / grid;
					
					// initial value
					double v0 [] = new double [] {5.0, 5.0, 5.0};
					
					// let go of the transient values
					for (int j = 0; j < 10000; ++j)
						v0 = rungeKutta (v0);
					
					xy [i] = getMaxima (v0);
				}
				repaint ();
			}
		});
	}

	// method to solve equation
	private double [] getMaxima (double [] v0)
	{
		// local variables
		double maxima [] = new double [iterate];
		double v1 [] = new double [v0.length];
		double v2 [] = new double [v0.length];
		
		for (int i = 0; i < iterate; ) {
			// evaluate new coordinate
			v1 = rungeKutta (v0);
			v2 = rungeKutta (v1);
			
			// check for fixed points
			if (Math.abs(v1[0] - v0[0]) < 0.0001 &&
					Math.abs(v1[0] - v2[0]) < 0.0001)
				return new double [] {v1[0]};
			
			// determine the z-maximum
			if (v1[0] > v0[0] && v1[0] > v2[0])
			{
				// assign maximum value
				maxima [i] = v1[0];
				
				// increment counter
				++i;
			}
			// set the new initial value
			v0 = v1;
		}
		
		return maxima;
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
		
		// van der pol equation
		vv[0] = -v[1] - v[2];
		vv[1] = v[0] + a * v[1];
		vv[2] = b + (v[0] - c) * v[2];
		
		return vv;
	}
	
	// method to paint on the panel
	public void paintComponent (Graphics g)
	{
		// check for multi ploting
		super.paintComponent(g);
			
		// make the draw area white
		g.setColor(Color.white);
		g.fillRect(xStart, yStart, width, height);
		g.setColor(Color.black);
		g.drawRect(xStart, yStart, width, height);
		
		// SHOW y-LEGEND
		g.drawString("" + yf, xStart - 25, (yStart + (int) (
				height * (yf - yf) / (yf - y0))));
		g.drawString("" + 0, xStart - 15, (yStart + (int) (
				height * (yf - 0) / (yf - y0))));
		g.drawString("" + y0, xStart - 25, (yStart + (int) (
				height * (yf - y0) / (yf - y0))));
		
		// SHOW x-legend
		g.drawString("" + x0, xStart - 15, yStart + height + 15);
		g.drawString("" + xf, xStart + width - 15, yStart + height + 15);
		
		for (int i = 0; i < 10; ++i)
			g.drawString("|", // + (x0 + i * (xf - x0) / 10), 
				xStart + i * width / 10, yStart + height + 5);
		
		// draw y-axis
		g.drawLine(
				(int) (xStart + -x0 / (xf - x0) * width), yStart + 10,
				(int) (xStart + -x0 / (xf - x0) * width), 
				yStart + height - 10);
		
		for (int i = 0; i < 10; ++i)
			g.drawString("-", xStart, yStart + i * height / 10 + 5);
		
		// draw x-axis
		g.drawString("" + xf, xStart + width - 15, yStart + height + 15);
		g.drawLine(xStart + 10, 
				(int) (yStart + yf / (yf - y0) * height),
				xStart + width - 10, 
				(int) (yStart + yf / (yf - y0) * height));
		
		// draw plot on the panel
		for (int i = 0; i < xy.length; ++i){
			// iterate through the vertical points
			for (int j = 0; j < xy[i].length; ++j) {
				g.drawOval ((int) (xStart + i * width / grid),
					(int)(yStart + (xy[i][j] - yf) / (y0 - yf) * height),
					1, 1);
			}
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		SwingUtilities.invokeLater (new Runnable () {
			public void run () {
				// instantiate the frame
				JFrame frame = new JFrame (
						"Bifurcation Diagram of Rossler Tent Map");
				frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
				frame.add(new RosslerBifurcation (), BorderLayout.CENTER);
				frame.setSize (screenWidth, screenHeight);
				frame.setLocationRelativeTo(null);
				frame.setVisible(true);
			}
		});
	}
}
