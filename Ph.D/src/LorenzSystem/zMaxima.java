package LorenzSystem;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Graphics;

import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;


public class zMaxima extends JPanel 
{
	// instance variables
	private static double y0 = 28.0;
	private static double yf = 50.0;
	private static double x0 = y0;
	private static double xf = yf;
	
	private int xStart = 50;
	private int yStart = 20;
	private int height = 500;
	private int width = 500;
	private static int screenWidth = 630;
	private static int screenHeight = 610;
	
	private double dt = 0.01;
	private double z[];
	private int n = 50000;
	
	// lorenz parameter
	private double rho = 28.0;
	private double sigma = 10.0;
	private double beta = 8.0 / 3.0;
	
	// Constructor
	public zMaxima()
	{
		z = new double [n];
		
		SwingUtilities.invokeLater(new Runnable () {
			public void run () {
				solveEquation (new double [] {5.0, 5.0, 5.0});
			}
		});
	}

	// method to solve equation
	private void solveEquation (double [] v0)
	{
		double v1 [] = new double [v0.length];
		double v2 [] = new double [v0.length];
		
		for (int i = 0; i < z.length; ) {
			// evaluate new coordinate
			v1 = rungeKutta (v0);
			v2 = rungeKutta (v1);
			
			// determine the z-maximum
			if (v1[2] > v0[2] && v1[2] > v2[2])
			{
				// assign maximum value
				z [i] = v1[2];
				
				// increment counter
				++i;
			}
			// set the new initial value
			v0 = v1;
		}
		// display the highes z-maxima
		double maximum = z[0];
		for (int i = 1; i < z.length; ++i)
			if (z[i] > maximum)
				maximum = z[i];
		
		JOptionPane.showMessageDialog(this, "z-maximum = " + maximum);
		
		repaint ();
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
		vv[0] = sigma * (v[1] - v[0]);
		vv[1] = v[0] * (rho - v[2]) - v[1];
		vv[2] = v[0] * v[1] - beta * v[2];
		
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
		g.drawString("" + 0, xStart + width / 2 - 5, yStart + height + 15);
		
		for (int i = 0; i < 10; ++i)
			g.drawString("" + (x0 + i), 
				xStart + i * width / 10 - 5, yStart + height + 15);
		
		// draw y-axis
		g.drawLine(
				(int) (xStart + -x0 / (xf - x0) * width), yStart + 10,
				(int) (xStart + -x0 / (xf - x0) * width), 
				yStart + height - 10);
		
		// draw x-axis
		g.drawString("" + xf, xStart + width - 15, yStart + height + 15);
		g.drawLine(xStart + 10, 
				(int) (yStart + yf / (yf - y0) * height),
				xStart + width - 10, 
				(int) (yStart + yf / (yf - y0) * height));
		
		// draw plot on the panel
		for (int i = 1; i < n; ++i){
			g.drawOval (
					(int) (xStart + (z[i - 1] - x0) / (xf - x0) * width),
					(int) (yStart + (z[i] - yf) / (y0 - yf) * height), 1, 1);
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
				JFrame frame = new JFrame ("Successive maxima of " +
						"z-coordinate of a Lorenz System");
				frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
				frame.add(new zMaxima (), BorderLayout.CENTER);
				frame.setSize (screenWidth, screenHeight);
				frame.setLocationRelativeTo(null);
				frame.setVisible(true);
			}
		});
	}
}
