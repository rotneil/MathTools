package nonlinearDyn;

import javax.swing.*;

import java.awt.*;
import java.awt.event.*;

public class KAMStandardMap extends JPanel 
	implements MouseListener, MouseMotionListener
{
	// instance variables
	private static double y0 = -Math.PI;
	private static double yf = Math.PI;
	/*
	private static double x0 = -3.0;
	private static double xf = 3.0;
	 */
	private static double x0 = -Math.PI;
	private static double xf = Math.PI;
	
	private int xStart = 50;
	private int yStart = 20;
	private int height = 450;
	private int width = 450;
	
	private double [] x;
	private double [] y;
	
	int n = 100000;
	private boolean multi_plot = false;
	
	private double a = 0.0;
	private double b = 1.0;
	
	// constructor
	public KAMStandardMap (boolean multi, 
			double x0, double xf, double y0, double yf)
	{
		// instantiate the instance variables
		multi_plot = multi;
		this.x0 = x0;
		this.y0 = y0;
		this.xf = xf;
		this.yf = yf;
		
		// prompt for the value of a
		try {
			a = Double.parseDouble(JOptionPane.showInputDialog(null, "Enter the a-value"));
			b = Double.parseDouble(JOptionPane.showInputDialog(null, "Enter the b-value"));
		} catch (NumberFormatException e) {}
		
		x = new double [n];
		y = new double [n];
		
		// register event handling
		addMouseListener (this);
		addMouseMotionListener (this);
		setBackground (java.awt.Color.white);
	}
	
	// method to solve equation
	private void solveEquation (double [] v0)
	{
		double v [] = new double [v0.length];
		v = v0;
		
		for (int i = 1; i < x.length; ++i) {
			v = S (v);
			
			x [i] = v [0];
			y [i] = v [1];
		}
		repaint ();
	}
	
	// equation definition
	private double [] S (double [] v)
	{
		double x, y;
		x = v[0] + v[1];
		y = b * v[1] + a * Math.sin(v[0] + v[1]);
		
		return new double [] {torus(x, 2.0 * Math.PI), torus(y, 2.0 * Math.PI)};
	}
	
	// bring x and y-values back to the region [-pi, pi]
	private double torus (double z, double mod)
	{
		int sign = (int) (z / Math.abs(z));
		int np = (int) (z / mod + sign * 0.5);
		z -= mod * np;
		
		return z;
	}
	
	// method to paint on the panel
	public void paintComponent (Graphics g)
	{
		// check for multi ploting
		if (!multi_plot) {
			super.paintComponent(g);
			
			// make the draw area white
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
		}
		
		// draw plot on the panel
		for (int i = 0; i < y.length ; ++i) {
			g.drawOval (
					(int) (xStart + (x[i] - x0) / (xf - x0) * width),
					(int) (yStart + (y[i] - yf) / (y0 - yf) * height), 1, 1);
		}
	}
	
	// method to launch the application
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		JFrame frame = new JFrame ("Standard Map Orbit at selected values of a");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(new KAMStandardMap (false, x0, xf, y0, yf));
		frame.setSize (580, 560);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
	}

	public void mouseDragged(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	public void mouseMoved(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	public void mouseClicked(MouseEvent e) {
		// TODO Auto-generated method stub
		// get click count
		switch (e.getClickCount()) {
		case 1:
			double xp = e.getX();
			double yp = e.getY();
			
			double xx = x0 + (xf - x0) * (xp - xStart) / width;
			double yy = yf + (y0 - yf) * (yp - yStart) / height;
			
			// initialize x- and y- coordinate
			y [0] = yy;
			x [0] = xx;
			
			multi_plot = true;
			solveEquation (new double [] {xx, yy});
			break;
		case 2:
			// prompt for the value of a
			try {
				a = Double.parseDouble(JOptionPane.showInputDialog(null, "Enter the a-value"));
				b = Double.parseDouble(JOptionPane.showInputDialog(null, "Enter the b-value"));
				
				multi_plot = false;
				solveEquation (new double [] {0, 0});
				
			} catch (NumberFormatException ee) {}
			
		}
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