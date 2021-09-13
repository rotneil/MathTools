package differentialEQ;

import javax.swing.*;

import java.awt.*;
import java.awt.event.*;

public class One_D_PhasePlot extends JPanel 
	implements MouseListener, MouseMotionListener
{
	// instance variables
	private static double a0 = 0.5;
	private static double y0 = -3.0;
	private static double yf = 3.0;
	private double x0;
	private double xf;
	
	
	private int xStart = 50;
	private int yStart = 20;
	private int height = 450;
	private int width = 450;
	
	private double [] x;
	private double [] y;
	private int grid = 2000;
	private double dt = 0.01;
	private boolean multi_plot = false;
	
	// constructor
	public One_D_PhasePlot (boolean multi, 
			double x0, double xf, double y0, double yf)
	{
		// instantiate the instance variables
		multi_plot = multi;
		this.x0 = x0;
		this.y0 = y0;
		this.xf = xf;
		this.yf = yf;
		x = new double [0];
		y = new double [0];
		
		// register event handling
		addMouseListener (this);
		addMouseMotionListener (this);
	}
	
	// method to solve equation
	private void solveEquation ()
	{
		for (int i = 1; i < x.length; ++i) {
			x [i] = x [i - 1] + dt;
			y [i] = rungeKutta (y[i - 1], x [i - 1]);
		}
		repaint ();
	}
	
	// Runge Kutta fourth order numerical integeration
	private double rungeKutta (double y, double t)
	{
		double c1, c2, c3, c4;
		
		c1 = f(y, t);
		
		c2 = y + dt * c1 / 2;
		c2 = f(c2, t + dt / 2);
		
		c3 = y + dt * c2 / 2;
		c3 = f(c3, t + dt / 2);
		
		c4 = y + dt * c3;
		c4 = f(c4, t + dt);
		
		c1 = y + dt * (c1 + 2 * (c2 + c3) + c4) / 6.0;
		
		return c1;
	}
	
	// equation definition
	private double f (double x, double t)
	{
		// logistic differential equation
		return (x * x * x - x);
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
		for (int i = 1; i < y.length && y [i] < yf && y [i] > y0; ++i) {
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
		JFrame frame = new JFrame ("The Logistic differential" +
				" model (a = ) " + a0 + " for different initial conditions");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(new One_D_PhasePlot (false, -3.0, 3.0, y0, yf));
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
		
		// convert xx to time factor
		int length = (int) ( (xf - xx) / dt);
		y = new double [length + 1];
		x = new double [length + 1];
		y [0] = yy;
		x [0] = xx;
		
		multi_plot = true;
		solveEquation ();
		
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