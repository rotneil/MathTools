package Pendulum;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.JFrame;
import javax.swing.JPanel;

public class Basin extends JPanel 
	implements MouseMotionListener, MouseListener {
	
	// instance variables
	// screen boundary values
	private int xStart = 55;
	private int yStart = 20;
	private int height = 450;
	private int width = 450;
	
	// plot dimension
	private double x0;
	private double xf;
	private double y0;
	private double yf;
	private int grid = 100;
	
	private int [] [] xy;
	
	// system variables
	private double rho = 1.66;
	private double c = 0.2;
	
	// diferential variables
	private static double dt = 0.01;
	private static double pi = Math.PI;
	private static double timeT = 2.0 * pi;
	
	// basin definition
	private int C = 0;
	private int B = 1;
	private int A = 2;
	private int NONE = 3;
	private double [] basinA = new double [] {-9.1291, -1.8577};
	private double [] basinB = new double [] {-0.7652, 0.3043};
	private double [] basinC = new double [] {-0.6565, -1.224};
	
	// variables used by mouse event
	private int mx1, mx2, my1, my2;
	private boolean draw = false;
	
	public Basin (final double x0, final double xf, 
			final double y0, final double yf)
	{
		this.x0 = x0;
		this.xf = xf;
		this.y0 = y0;
		this.yf = yf;
		
		xy = new int [grid] [grid];
		
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {

			public void run() 
			{
				// local variables of the grid point
				double [] v = new double [2];
				
				for (int i = 0; i < grid; ++i) 
				{
					// assign the horizontal value
					v[0] = x0 + i * (xf - x0) / grid;
					
					// assign the vertical value
					for (int j = 0; j < grid; ++j) 
					{
						v[1] = y0 + j * (yf - y0) / grid;
						
						// perform 500-iterates of the point
						for (int k = 0; k < 1000; ++k)
							v = fT(v);
						
						// decide to which basin v has been trapped
						xy[i] [j] = findBasin (v);
					}
					// show the basin
					repaint ();
				}	// end for loop
				
				repaint ();
			}	// end method run
			
		});	// end executor
		executor.shutdown();
		
		// add mouse motion listener
		setBackground (Color.white);
		addMouseListener (this);
		addMouseMotionListener (this);
	}
	
	// method to declare the basin to which an iterate has been trapped
	private int findBasin (double [] v)
	{
		if (Math.abs(v[0] - basinA[0]) <= 0.0001 && Math.abs(v[1] - basinA[1]) <= 0.0001)
			return A;
		else if (Math.abs(v[0] - basinB[0]) <= 0.0001 && Math.abs(v[1] - basinB[1]) <= 0.0001)
			return B;
		else if (Math.abs(v[0] - basinC[0]) <= 0.0001 && Math.abs(v[1] - basinC[1]) <= 0.0001)
			return C;
		else
			return NONE;
	}
	
	// define a timeT function of the system
	private double [] fT (double [] v)
	{
		// define the start time
		double [] vv = new double [v.length];
		vv = v;
		
		// iterate the system for timeT
		for (double t = 0; t < timeT; t += dt) {
			vv = rungeKutta (vv, t);
			
			// bring back the value of x
			int sign = (int) (vv[0] / Math.abs(vv[0]));
			int np = (int) (vv[0] / (2.0 * pi) + sign * 0.5);
			vv[0] -= np * 2.0 * pi;
		}
		
		return vv;
	}
	
	// method rungeKutta
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
	
	// define the system equation
	private double [] f (double [] v, double t)
	{
		// local variables
		double [] vv = new double [v.length];
		
		// assign values of vv
		vv [0] = v [1];
		vv [1] = rho * Math.sin(t) - c * v [1] - Math.sin(v[0]);
		
		return vv;
	}
	
	// method to paint on the panel
	public void paintComponent (Graphics g)
	{
		super.paintComponent(g);
		
		// draw plot on the panel
		for (int i = 0; i < xy.length; ++i) {
			int xx = xStart + (int) (i * width / grid);
			for (int j = 0; j < xy [i].length; ++j) {
				int yy = yStart + height - (int) (j * height / grid);
				
				// check for the basin
				switch (xy[i][j]) {
					case 0:	// A
						g.setColor(Color.white);
						g.drawOval(xx, yy, 1, 1);
						break;
					case 1:	// B
						g.setColor(Color.black);
						g.drawOval(xx, yy, 1, 1);
						break;
					case 2:	// C
						g.setColor(Color.gray);
						g.drawOval(xx, yy, 1, 1);
						break;
					default:
						g.setColor(Color.red);
						g.drawOval(xx, yy, 1, 1);
				}
				
			}
		}
		
		// draw the selection area
		g.setColor(Color.black);
		if (draw)
			g.drawRect(Math.min(mx1, mx2), Math.min(my1, my2),
					Math.abs(mx1 - mx2), Math.abs(my1 - my2));
		
		// make the draw area white
		g.drawRect(xStart, yStart, width, height);
		
		// show the legend
		g.drawString("" + todp(yf, 4), xStart - 25, yStart + 5);
		g.drawString("" + todp(y0, 4), xStart - 25, yStart + height + 5);
		g.drawString("" + todp(x0 / pi, 4) + "pi", xStart - 5, yStart + height + 15);
		g.drawString("" + todp(xf / pi, 4) + "pi", xStart + width - 10, yStart + height + 15);
		
		// show minor marks
		for (int i = 0; i < 21; ++i) {
			// divide y-axis
			g.drawString("-", xStart, yStart + (int) (i * (height) / 20) + 4);
			g.drawString("-", xStart + width - 2, yStart + (int) (i * (height) / 20) + 4);
			g.drawString("'", xStart + (int) (i * width / 20) - 1, yStart + 10);
			g.drawString("'", xStart + (int) (i * width / 20) - 1, yStart + (height) + 7);
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		double x0, xf, y0, yf;
		x0 = -pi;
		xf = pi;
		y0 = -2.0;
		yf = 4.0;
		// TODO Auto-generated method stub
		JFrame frame = new JFrame ("The Basin of forced damped pendulum for [-pi, pi] X ["
				+ y0 + ", " + yf + "]");
		//frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(new Basin (x0, xf, y0, yf));
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
	private double todp (double x, int dp)
	{
		double sf = Math.floor(Math.pow(10, dp) * x + 0.5);
		return sf / Math.pow(10, dp);
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
		x0 = todp (x0, 4);
		xf = todp (xf, 4);
		y0 = todp (y0, 4);
		yf = todp (yf, 4);
		
		// TODO Auto-generated method stub
		JFrame frame = new JFrame ("The bifurcation Diagram of "
				+ "logistic map for [" + todp(x0 / pi, 4) + "pi, " + todp(xf / pi, 4) + "pi] X ["
				+ y0 + ", " + yf + "]");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(new Basin (x0, xf, y0, yf));
		frame.setSize (580, 560);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
	}
}