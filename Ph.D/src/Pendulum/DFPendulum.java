package Pendulum;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

/*
** A program to study the driven pendulum under damping 
** via the fourth-order Runge-Kutta algorithm
**
*/

public class DFPendulum extends JPanel implements MouseListener, MouseMotionListener, KeyListener
{
	// instance variables
	private static double y0 = -10.0;
	private static double yf = 10.0;
	private static double x0 = -Math.PI;
	private static double xf = Math.PI;
	
	private int xStart = 50;
	private int yStart = 20;
	private int height = 450;
	private int width = 450;
	
	private double [] x;
	private double [] y;
	
	private double dt = 0.01;
	double c = 0.2, rho = 1.66;
	
	private boolean multi_plot = false;
	private int keyN = 0;
	private double steadyX, steadyY;
	
	// constructor
	public DFPendulum (boolean multi, double x0, double xf, double y0, double yf)
	{
		// instantiate the instance variables
		multi_plot = multi;
		this.x0 = x0;
		this.y0 = y0;
		this.xf = xf;
		this.yf = yf;
		x = new double [5000];
		y = new double [5000];
		
		// register event handling
		addMouseListener (this);
		addMouseMotionListener (this);
	}
	
	// method to solve equation
	private void solveEquation (double [] v0)
	{
		double [] v = new double [v0.length];
		v = v0;
		
		for (int i = 0; i < x.length; ++i) {
			v = rungeKutta (v, i * dt);
			v [0] = torus (v[0], 2.0 * Math.PI);
			
			x [i] = v [0];
			y [i] = v [1];
		}
		repaint ();
	}
	
	// equation definition
	private double [] f (double [] v, double t)
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
	public void paintComponent (Graphics g)
	{
		// check for multi ploting
		if (!multi_plot) {
			super.paintComponent(g);
			
			// make the draw area white
			setBackground(java.awt.Color.white);
			g.drawRect(xStart, yStart, width, height);
			
			// SHOW y-LEGEND
			g.drawString("" + yf, xStart - 25, (yStart + (int) (
					height * (yf - yf) / (yf - y0))));
			g.drawString("" + 0, xStart - 15, (yStart + (int) (
					height * (yf - 0) / (yf - y0))));
			g.drawString("" + y0, xStart - 25, (yStart + (int) (
					height * (yf - y0) / (yf - y0))));
			
			// SHOW x-legend
			g.drawString(todp (x0 / Math.PI, 0) + "pi", xStart - 15, yStart + height + 15);
			g.drawString(todp (xf / Math.PI, 0) + "pi", xStart + width - 15, yStart + height + 15);
			g.drawString("" + 0, xStart + width / 2 - 5, yStart + height + 15);
			
			// show the x and y-axes
			g.drawLine(xStart + width / 2, yStart + 10,
					xStart + width / 2, yStart + height - 10);
			g.drawLine(xStart + 10, yStart + height / 2,
					xStart + width - 10, yStart + height / 2);
		}
		
		// draw plot on the panel
		for (int i = 0; i < x.length; ++i)
			g.drawOval (
					(int) (xStart + (x[i] - x0) / (xf - x0) * width),
					(int) (yStart + (y[i] - yf) / (y0 - yf) * height), 1, 1);
		/*
			g.drawLine (
					(int) (xStart + (x[i - 1] - x0) / (xf - x0) * width),
					(int) (yStart + (y[i - 1] - yf) / (y0 - yf) * height),
					(int) (xStart + (x[i] - x0) / (xf -x0) * width),
					(int) (yStart + (y[i] - yf) / (y0 - yf) * height));
		*/
	}
	
	// method to convert values to certain decimal places
	private double todp (double value, int p)
	{
		// get the sign of value
		int sign = value > 0 ? 1 : -1;
		int vv = (int) (value * Math.pow(10, p) + sign * 0.5);
		return vv / Math.pow(10, p);
	}
	
	// method to launch the application
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		JFrame frame = new JFrame ("Forced Damped Pendulum with different initial condition");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		DFPendulum pendulum = new DFPendulum (false, x0, xf, y0, yf);
		frame.add(pendulum);
		frame.addKeyListener(pendulum);
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
		x = new double [5000];
		y = new double [5000];
		
		// use click count to 
		switch (e.getClickCount()) {
			case 1:
				if (xx >= x0 && xx <= xf && yy >= y0 && yy <= yf) {
					multi_plot = true;
					solveEquation (new double [] {xx, yy});
				}
				break;
			case 2:
				// prompt for rho and c
				try {
					rho = Double.parseDouble(
							JOptionPane.showInputDialog(null, "Enter the value of rho"));
					c = Double.parseDouble(
							JOptionPane.showInputDialog(null, "Enter the value of c"));
					
					// change plot param
					multi_plot = false;
					solveEquation (new double [] {xx, yy});
				} catch (NumberFormatException ee){} catch (NullPointerException ee) {}
			}
		steadyX = x [4999];
		steadyY = y [4999];
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

	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub
		// Key_C for clear
		// key_o for one increment
		// Key_T for ten increment
		// key_H for hundred
		// key_n for thousand increment and 
		// key_m for ten thousand
		switch (e.getKeyChar()) 
		{
			case 'c' : case 'C':
				keyN = 1;
				break;
			case 'o' : case 'O':
				++keyN;
				break;
			case 't' : case 'T':
				keyN += 10;
				break;
			case 'h' : case 'H':
				keyN += 100;
				break;
			case 'n' : case 'N':
				keyN += 1000;
				break;
			case 'm' : case 'M':
				keyN += 10000;
				break;
			default:
		}
		// change plot param
		x = new double [keyN];
		y = new double [keyN];
		
		multi_plot = false;
		solveEquation (new double [] {steadyX, steadyY});
	}

	public void keyPressed(KeyEvent e) {
		// TODO Auto-generated method stub
		
	}

	public void keyReleased(KeyEvent e) {
		// TODO Auto-generated method stub
		
	}
}