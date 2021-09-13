package nonlinearDyn;

import java.awt.Color;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class IteratedFunctionalSystem extends javax.swing.JPanel 
{
	private static final long serialVersionUID = 1L;

	// instance variables
	private static int width = 500;
	private static int height = 500;
	private static int xStart = 50;
	private static int yStart = 25;
	
	private double x0, xf, y0, yf;
	private int n = 1000000;
	private double [] [] ifs;
	private double [] x, y;
	private java.util.Random random;
	
	// no-argument constructor
	public IteratedFunctionalSystem (double x0, double xf, double y0, double yf) 
	{
		this.x0 = x0;
		this.xf = xf;
		this.y0 = y0;
		this.yf = yf;
		
		// instantiate the iterated functional system for barnsley fern
		x = new double [n];
		y = new double [n];
		
		ifs = new double [] [] {
				{0.0, 0.0, 0.0, 0.16, 0.0, 0.0, 0.01},
				{0.85, 0.04, -0.04, 0.85, 0.0, 1.60, 0.85},
				{0.20, -0.26, 0.23, 0.22, 0.0, 1.60, 0.07},
				{-0.15, 0.28, 0.26, 0.24, 0.0, 0.44, 0.07}
		};
		/*
		ifs = new double [] [] {
				{0.0, 0.0, 0.0, 0.25, 0.0, -0.4, 0.02},
				{0.95, 0.005, -0.005, 0.93, -0.002, 0.5, 0.84},
				{0.035, -0.2, 0.16, 0.04, -0.09, 0.02, 0.07},
				{-0.04, 0.2, 0.16, 0.04, 0.083, 0.12, 0.07}
		};*/
		
		random = new java.util.Random();

		x[0] = 0.0;
		y[0] = 0.0;
		
		// called the fractal generator
		generateFractal ();
	}
	
	// method to iteratively generate the fractal
	private void generateFractal ()
	{
		// execute in a different thread
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {
			@Override
			public void run () {
				// define the ifs probability
				int [] p = new int [ifs.length];
				p[0] = (int) (ifs[0][6] * 100);
				
				for (int i = 1; i < p.length; ++i)
					p[i] = p[i - 1] + (int) (ifs[i][6] * 100);
				
				// randomly populate the x- and y-cordinates
				int level = 1;
				while (level < n) {
					int ran = random.nextInt(100);
					
					// choose the ifs to use
					int f = 0;
					for (int i = 0; i < ifs.length; ++i) {
						if (p[i] > ran) {
							f = i;
							break;
						}
					}
					
					// manipulate (x, y) based on the chosen ifs
					x[level] = ifs [f] [0] * x [level - 1] +
							ifs [f] [1] * y [level - 1] + ifs [f][4];
					y[level] = ifs [f] [2] * x [level - 1] +
							ifs [f] [3] * y [level - 1] + ifs [f][5];
					
					// increment level
					++level;
				}	// end while block
				
				// repaint
				javax.swing.SwingUtilities.invokeLater(new Runnable () {
					@Override
					public void run () {
						repaint ();
					}
				});
			}
		});
	}	// end method generate Fractal
	
	
	// method to paint fractal on the screen
	@Override
	public void paintComponent (java.awt.Graphics g)
	{
		super.paintComponent (g);
		setBackground (java.awt.Color.white);
		
		// draw the frame
		//g.setColor(java.awt.Color.blue);
		//g.drawRect(xStart, yStart, width, height);
		
		// set the paint color to blue and start to draw
		g.setColor(new java.awt.Color (55, 159, 85));
		for (int i = 0; i < n; ++i)
			g.drawOval(
					(int) (xStart + (x[i] - x0) / (xf - x0) * width),
					(int) (yStart + (y[i] - yf) / (y0 - yf) * height), 1, 1);
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		// instantiate the window frame
		javax.swing.JFrame frame = new javax.swing.JFrame(
				"Intrated Functional System");
		frame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
		frame.add(new IteratedFunctionalSystem (-3.0, 3.0, 0.0, 10.0));
		frame.setSize(630, 610);
		frame.setVisible(true);
	}

}
