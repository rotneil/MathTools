package duffingOscillator;

import java.awt.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.*;

public class DuffingOscillatorBasin extends JPanel
{
	// instance variables
	private static double y0 = -10.0;
	private static double yf = 10.0;
	private static double x0 = -10.0;
	private static double xf = 10.0;
	
	private int xStart = 50;
	private int yStart = 20;
	private int height = 500;
	private int width = 500;
	private static int screenWidth = 630;
	private static int screenHeight = 610;
	
	private int grid = 500;
	private double dt = 0.01;
	private int iterate = 40000;
	private int [] [] xy;
	
	private static int INFINITY = 0;
	private static int BASIN_A = 1;
	private static int BASIN_B = 2;
	
	// constructor
	public DuffingOscillatorBasin ()
	{
		xy = new int [grid + 1] [grid + 1];
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {
			public void run () {
				drawBasin ();
			}
		});
	}
	
	// method to populate the grids
	public void drawBasin ()
	{
		// iterate through the initial conditions
		for (int i = 0; i < xy.length; ++i) {
			double xx = x0 + i * (xf - x0) / grid;
			
			for (int j = 0; j < xy [i].length; ++j) {
				double yy = y0 + j * (yf - y0) / grid;
				
				double [] v = new double [] {xx, yy};
				
				// solve each initial condition
				for (int k = 0; k < iterate; ++k)
					v = rungeKutta (v, v[0]);
				
				if (Math.abs(v[1]) < 0.0001) {
					if (Math.abs(v[0] + 1.0) < 0.0001)
						xy [i] [j] = BASIN_A;
					else if (Math.abs(v[0] - 1.0) < 0.0001)
						xy [i] [j] = BASIN_B;
				}
			}
		}
		repaint ();
	}
	
	// method to paint component
	public void paintComponent (Graphics g)
	{
		// check if to draw the basin
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
		
		// draw the basin
		for (int i = 0; i < xy.length; ++i) {
			int xx = (int) (xStart + i * width / grid);
			
			for (int j = 0; j < xy.length; ++j) {
				int yy = (int) (yStart + height * (grid - j) / grid);
				
				// choose basin color
				if (xy [i] [j] == BASIN_A) {
					g.setColor(Color.black);
					g.drawOval(xx, yy, 1, 1);
				} else if (xy [i] [j] == BASIN_B) {
					g.setColor(Color.lightGray);
					g.drawOval(xx, yy, 1, 1);
				}
			}
		}
	}
	
	// method to launch the application
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		JFrame frame = new JFrame (
				"Basin of Forced-damped Duffing Oscillator " +
					"[" + x0 + ", " + xf + "] X [" + y0 + ", " + yf + "]");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(new DuffingOscillatorBasin ());
		frame.setSize (screenWidth, screenHeight);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
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
	
	// equation definition
	private double [] f (double [] v, double t)
	{
		int l = v.length;
		double [] vv = new double [l];
		
		// duffing equation
		vv[0] = v[1];
		vv[1] = v[0] - Math.pow(v[0], 3) - 0.1 * v[1];
		
		return vv;
	}
}