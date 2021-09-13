package henonMap;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.JFrame;
import javax.swing.JPanel;

public class Bifurcation extends JPanel 
	implements MouseMotionListener, MouseListener {
	
	// instance variables
	private int xStart = 50;
	private int yStart = 20;
	private int height = 450;
	private int width = 450;
	
	private double x0;
	private double xf;
	private double y0;
	private double yf;
	private int grid = 1000;
	private int yGrid = 500;
	
	private double [] [] xy;
	private double a;
	private static double b = 0.3;
	
	// variables used by mouse event
	private int mx1, mx2, my1, my2;
	private boolean draw = false;
	
	public Bifurcation (final double x0, final double xf, 
			final double y0, final double yf)
	{
		this.x0 = x0;
		this.xf = xf;
		this.y0 = y0;
		this.yf = yf;
		
		xy = new double [grid] [yGrid];
		
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {

			public void run() 
			{
				// local point definition
				double [] v = new double [2];
				
				for (int i = 0; i < grid; ++i) 
				{
					a = x0 + i * (xf - x0) / grid;
					
					// choose random value for x and y
					v [0] = v [1] = 0.5;
					
					// perferm k-iterates to move the point to stability (if there exists)
					for (int k = 0; k < 100; ++k)
						v = f(v);
					
					// perform bifurcation plot
					for (int j = 0; j < xy[i].length; ++j) {
						v = f (v);
						xy [i] [j] = v[0];
					}
				}	// end for loop
				
				repaint ();
			}	// end method run
			
		});	// end executor
		
		// add mouse motion listener
		addMouseListener (this);
		addMouseMotionListener (this);
	}
	
	// define the system equation
	private double [] f (double [] v)
	{
		// local variables
		double [] vv = new double [v.length];
		
		vv [0] = a - sq(v[0]) + b * v[1];
		vv [1] = v[0];
		
		// define variables;
		return vv;
	}
	
	private double sq (double x) {return x * x;}
	
	// method to paint on the panel
	public void paintComponent (Graphics g)
	{
		super.paintComponent(g);
		
		// make the draw area white
		g.setColor(Color.white);
		g.fillRect(xStart, yStart, width, height);
		g.setColor(Color.black);
		g.drawRect(xStart, yStart, width, height);
		
		// show the legend
		g.setColor(java.awt.Color.black);
		g.drawString("" + yf, xStart - 25, yStart + 5);
		g.drawString("" + y0, xStart - 25, yStart + height + 5);
		g.drawString("" + x0, xStart - 5, yStart + height + 15);
		g.drawString("" + xf, xStart + width - 10, yStart + height + 15);
		
		// show minor marks
		for (int i = 0; i < 21; ++i) {
			// divide y-axis
			g.drawString("-", xStart, yStart + (int) (i * height / 20) + 5);
			g.drawString("-", xStart + width - 2, yStart + (int) (i * height / 20) + 5);
			g.drawString("'", xStart + (int) (i * width / 20), yStart + 10);
			g.drawString("'", xStart + (int) (i * width / 20), yStart + height + 7);
		}
		
		// draw plot on the panel
		for (int i = 0; i < xy.length; ++i) {
			int xx = xStart + (int) (i * width / xy.length);
			for (int j = 0; j < xy [i].length; ++j) {
				int yy = yStart + (int) (
						height * (yf - xy [i] [j]) / (yf - y0) );
				
				// draw a point here
				g.drawOval(xx, yy, 1, 1);
			}
		}
		
		if (draw)
			g.drawRect(Math.min(mx1, mx2), Math.min(my1, my2),
					Math.abs(mx1 - mx2), Math.abs(my1 - my2));
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		double x0, xf, y0, yf;
		x0 = -1.0;
		xf = 2.0;
		y0 = -3.0;
		yf = 4.0;
		// TODO Auto-generated method stub
		JFrame frame = new JFrame ("The bifurcation Diagram of "
				+ "Henon map for [" + x0 + ", " + xf + "] X ["
				+ y0 + ", " + yf + "] for b = " + b);
		//frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(new Bifurcation (x0, xf, y0, yf));
		frame.setSize (580, 560);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
	}

	public void mouseClicked(MouseEvent e) {
		// TODO Auto-generated method stub
		draw = false;
	}

	public void mousePressed(MouseEvent e) {
		// TODO Auto-generated method stub
		mx1 = e.getX();
		my1 = e.getY();
		draw = false;
	}

	// method to give a 2 s.f double
	private double to4sf (double x)
	{
		double sf = Math.floor(10000 * x + 0.5);
		return sf / 10000.0;
	}
	
	public void mouseReleased(MouseEvent e) {
		// TODO Auto-generated method stub
		// start another frame
		if (draw)
			startAnotherFrame (
				getXValue (Math.min(mx1, mx2)),
				getXValue (Math.max(mx1, mx2)),
				getYValue (Math.max(my1, my2)),
				getYValue (Math.min(my1, my2)));
	}

	public void mouseEntered(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	public void mouseExited(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	public void mouseDragged(MouseEvent e) {
		// TODO Auto-generated method stub
		mx2 = e.getX();
		my2 = e.getY();
		draw = true;
	}

	public void mouseMoved(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	private double getXValue (int xx)
	{
		return (x0 + (xf - x0) * (xx - xStart) / width);
	}
	
	private double getYValue (int yy)
	{
		return (yf + (y0 - yf) * (yy - yStart) / height);
	}
	
	private void startAnotherFrame (
			double x0, double xf, double y0, double yf)
	{
		x0 = to4sf(x0);
		xf = to4sf(xf);
		y0 = to4sf(y0);
		yf = to4sf(yf);
		
		// TODO Auto-generated method stub
		JFrame frame = new JFrame ("The bifurcation Diagram of "
				+ "Henon Map [" + x0 + ", " + xf + "] X ["
				+ y0 + ", " + yf + "] for a = " + a + " and b = " + b);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(new Bifurcation (x0, xf, y0, yf));
		frame.setSize (580, 560);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
	}
}