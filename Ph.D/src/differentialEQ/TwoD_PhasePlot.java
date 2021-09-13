package differentialEQ;

import javax.swing.*;

import java.awt.*;
import java.awt.event.*;

public class TwoD_PhasePlot extends JPanel 
	implements MouseListener, MouseMotionListener
{
	// instance variables
	private static double y0 = -10.0;
	private static double yf = 10.0;
	private static double x0 = -3.0;
	private static double xf = 3.0;
	
	private int xStart = 50;
	private int yStart = 20;
	private int height = 450;
	private int width = 450;
	
	private double [] t;
	private double [] x;
	private double [] y;
	
	int n = 10000;
	private double dt = 0.01;
	private boolean multi_plot = false;
	
	// constructor
	public TwoD_PhasePlot (boolean multi, 
			double x0, double xf, double y0, double yf)
	{
		// instantiate the instance variables
		multi_plot = multi;
		this.x0 = x0;
		this.y0 = y0;
		this.xf = xf;
		this.yf = yf;
		t = new double [n];
		x = new double [n];
		y = new double [n];
		
		// register event handling
		addMouseListener (this);
		addMouseMotionListener (this);
	}
	
	// method to solve equation
	private void solveEquation (double [] v0)
	{
		double v [] = new double [v0.length];
		v = v0;
		
		for (int i = 1; i < x.length; ++i) {
			t [i] = t [i - 1] + dt;
			v = rungeKutta (v, x [i - 1]);
			
			x [i] = v [0];
			y [i] = v [1];
		}
		repaint ();
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
		
		/*// real and unique eigenvalues
		vv [0] = 1 * v[0] - 2 * v[1];
		vv [1] = 5 * v[0] - 1 * v[1];
		*/
		
		/*// real and repeated eigenvalues
		vv[0] = 3 * v[0] + v[1];
		vv[1] = 3 * v[1];
		*/
		
		/*//imaginary eigenvalues
		vv[0] = v[1];
		vv[1] = -v[0];
		*/
		
		/*// complex roots
		double a = 1.0;
		double b = 10.0;
		vv[0] = a * v[0] - b * v[1];
		vv[1] = b * v[0] + a * v[1];
		*/
		
		// nonlinear system
		vv[0] = v[1] - v[0] * v[0];
		vv[1] = -v[0];
		
		return vv;
	}
	
	// method to paint on the panel
	public void paintComponent (Graphics g)
	{
		// check for multi ploting
		if (!multi_plot) {
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
			
			// show the x and y-axes
			g.drawLine(xStart + width / 2, yStart + 10,
					xStart + width / 2, yStart + height - 10);
			g.drawLine(xStart + 10, yStart + height / 2,
					xStart + width - 10, yStart + height / 2);
		}
		
		// draw plot on the panel
		for (int i = 1; i < y.length && y [i] < yf && y [i] > y0
				&& x [i] >= x0 && x [i] <= xf; ++i) {
			g.drawLine (
					(int) (xStart + (x[i - 1] - x0) / (xf - x0) * width),
					(int) (yStart + (y[i - 1] - yf) / (y0 - yf) * height),
					(int) (xStart + (x[i] - x0) / (xf -x0) * width),
					(int) (yStart + (y[i] - yf) / (y0 - yf) * height));
		}
	}
	
	// method to launch the application
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		JFrame frame = new JFrame ("Two D phase Plot for different initial conditions");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(new TwoD_PhasePlot (false, x0, xf, y0, yf));
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
		double xp = e.getX();
		double yp = e.getY();
		
		double xx = x0 + (xf - x0) * (xp - xStart) / width;
		double yy = yf + (y0 - yf) * (yp - yStart) / height;
		
		// initialize x- and y- coordinate
		y [0] = yy;
		x [0] = xx;
		
		multi_plot = true;
		solveEquation (new double [] {xx, yy});
		
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