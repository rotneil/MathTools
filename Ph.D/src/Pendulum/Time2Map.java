package Pendulum;

import javax.swing.*;

import java.awt.*;
import java.awt.event.*;

public class Time2Map extends JPanel 
	implements MouseListener, MouseMotionListener
{
	// instance variables
	private static double y0 = -6.5;
	private static double yf = 1.5;
	private static double x0 = -Math.PI;
	private static double xf = Math.PI;
	
	private int xStart = 50;
	private int yStart = 20;
	private int height = 500;
	private int width = 500;
	private static int screenWidth = 630;
	private static int screenHeight = 610;
	
	private double [] x;
	private double [] y;
	
	int n = 1000000;
	private double dt = 0.01;
	private int time2pi = (int) (2.0 * Math.PI / dt);
	private static double rho = 2.5;
	private static double c = 0.05;
	
	// constructor
	public Time2Map ()
	{
		// instantiate the instance variables
		x = new double [n];
		y = new double [n];
		
		// register event handling
		addMouseListener (this);
		addMouseMotionListener (this);
	}
	
	// method to solve equation
	private void solveEquation (final double [] v0)
	{
		SwingUtilities.invokeLater (new Runnable () {
			public void run () {
				double v [] = new double [v0.length];
				v = v0;
				
				// discard the first 1000 orbits
				for (int i = 0; i < 1000; ++i)
					v = rungeKutta (v);
				
				// iterate to make n-number of points for the attractor
				for (int i = 0; i < x.length; ++i) {
					// get point each time2pi
					for (int k = 0; k < time2pi; ++k) {
						v = rungeKutta (v);
						
						// restrict theta to with +pi and -pi
						int sign = v[0] >= 0 ? 1 : -1;
						int np = (int) (v[0] / (2 * Math.PI) + sign * 0.5);
						v[0] -= 2 * np * Math.PI;
					}
					x[i] = v[0];
					y[i] = v[1];
					v[2] = 0.0; 	// implement the time-2pi
				}
				
				repaint ();
			}
		});
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
		double [] vv = new double [v.length];
		
		// duffing equation
		vv[0] = v[1];
		vv[1] = -c * v[1] - Math.sin(v[0]) + rho * Math.sin(v[2]);
		vv[2] = 1.0;
		
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
		
		// draw plot on the panel
		for (int i = 0; i < y.length; ++i){
			g.drawOval (
				(int) (xStart + (x[i] - x0) / (xf - x0) * width),
				(int) (yStart + (y[i] - yf) / (y0 - yf) * height), 1, 1);
		}
	}
	
	// method to launch the application
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		// launch application
		SwingUtilities.invokeLater(new Runnable () {
			public void run () {
				// declare the frame to house the panel
				JFrame frame = new JFrame ();
				
				// prompt for the value of c
				try {
					c = Double.parseDouble(JOptionPane.showInputDialog(
							frame, "Enter the value of c"));
				} catch (NumberFormatException e) {}
				
				frame.setTitle("Time2pi Map of Forced Damped " +
						"Pendulum with c = " +	c + " and rho = " + rho);
				frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
				frame.add(new Time2Map (),
						BorderLayout.CENTER);
				frame.setSize (screenWidth, screenHeight);
				frame.setLocationRelativeTo(null);
				frame.setVisible(true);
			}
		});
	}

	public void mouseDragged(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	public void mouseMoved(MouseEvent e) {
		// TODO Auto-generated method stub
	}
	
	public void mouseClicked(MouseEvent e) {
		// TODO Auto-generated method stub
		double xp = e.getX();
		double yp = e.getY();
		
		double xx = x0 + (xf - x0) * (xp - xStart) / width;
		double yy = yf + (y0 - yf) * (yp - yStart) / height;
		
		// initialize x- and y- coordinate
		y [0] = yy;
		x [0] = xx;
		
		solveEquation (new double [] {xx, yy, 0.0});
		
	}

	public void mousePressed(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	public void mouseReleased(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	public void mouseEntered(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	public void mouseExited(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
	
}