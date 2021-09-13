package york.challenge7;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;
import javax.swing.Timer;


public class Step2 extends JPanel 
implements MouseListener, MouseMotionListener
{
	// instance variables
	private static double y0 = -4.0;
	private static double yf = 4.0;
	private static double x0 = -4.0;
	private static double xf = 4.0;
	
	private int xStart = 50;
	private int yStart = 20;
	private int height = 500;
	private int width = 500;
	private static int screenWidth = 630;
	private static int screenHeight = 610;
	
	private double dt = 0.01;
	private double x[], y[];
	private int n = (int) (2 * Math.PI / dt);
	private int nStep;
	
	// Constructor
	public Step2()
	{
		x = new double [n];
		y = new double [n];
		
		addMouseMotionListener (this);
		addMouseListener (this);
		
		ActionListener animator = new ActionListener () {

			public void actionPerformed(ActionEvent e) {
				// TODO Auto-generated method stub
				for (nStep = 0; nStep < n; nStep += 20)
					repaint ();
			}
		};
		new Timer (1000, animator).start();
	}

	// method to solve equation
	private void solveEquation (double [] v0)
	{
		double v [] = new double [v0.length];
		v = v0;
		
		for (int i = 1; i < x.length; ++i) {
			v = rungeKutta (v);
			
			x [i] = v [0];
			y [i] = v [1];
		}
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
		vv[0] = v[1] - (Math.pow(v[0], 3) - v[0]);
		vv[1] = -v[0];
		
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
		
		// plot the function F(x)
		for (double x = x0; x < xf; x += 0.05) {
			double xxx = x + 0.05;
			g.drawLine(
					(int) (xStart + (x - x0) / (xf - x0) * width),
					(int) (yStart + (F(x) - yf) / (y0 - yf) * height),
					(int) (xStart + (xxx - x0) / (xf - x0) * width),
					(int) (yStart + (F(xxx) - yf) / (y0 - yf) * height));
		}
		
		// draw plot on the panel
		for (int i = 1; i < nStep && i < n; ++i){
			if (y [i] <= yf && y [i] >= y0 && x [i] >= x0 
					&& x [i] <= xf)
				g.drawLine (
					(int) (xStart + (x[i - 1] - x0) / (xf - x0) * width),
					(int) (yStart + (y[i - 1] - yf) / (y0 - yf) * height),
					(int) (xStart + (x[i] - x0) / (xf -x0) * width),
					(int) (yStart + (y[i] - yf) / (y0 - yf) * height));
		}
	}
	
	private double F (double x)
	{
		return (Math.pow(x, 3) / 3.0 - x);
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		SwingUtilities.invokeLater (new Runnable () {
			public void run () {
				// instantiate the frame
				JFrame frame = new JFrame ("Van der Pol Challenge");
				frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
				frame.add(new Step2 (), BorderLayout.CENTER);
				frame.setSize (screenWidth, screenHeight);
				frame.setLocationRelativeTo(null);
				frame.setVisible(true);
			}
		});
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
		
		solveEquation (new double [] {xx, yy});
	}
	
	public void mouseDragged(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	public void mouseMoved(MouseEvent e) {
		// TODO Auto-generated method stub
		
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
