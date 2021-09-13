/*
** This program plots the stable and unstable manifold of 
** Henon Map defined by
** 		f(x, y) = (a - x*x + by, x)
**
*/

package manifold;

import java.awt.Color;
import java.awt.Graphics;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;


public class HenonMap extends JPanel
{
	// instance variables
	private int xStart = 50;
	private int yStart = 20;
	private int width = 500;
	private int height = 500;
	
	private double x0;
	private double xf;
	private double y0;
	private double yf;
	
	// henon variables
	private double a = 2.12;
	private double b = -0.3;
	
	private double epsilon = 0.000001;
	
	private double [] [] [] stable = new double [1000] [1000] [2];
	private double [] [] unstable = new double [1000] [2];
	
	// constructor
	public HenonMap (double x0, double xf, double y0, double yf)
	{
		// set the instance variables
		this.x0 = x0;
		this.xf = xf;
		this.y0 = y0;
		this.yf = yf;
		
		// evaluate the saddle fixed point and eigenvector
		double [] p = new double [] {0.944522, 0.944522};
		double [] vU = new double [] {1.0, -0.58};
		double [] vS = new double [] {1.0, -5.71339};
		
		// evaluate the 10^-6 of p
		double [] ppU = new double [] {p[0] + epsilon * vU[0], p[1] + epsilon * vU[1]};
		double [] ppS = new double [] {p[0] + epsilon * vS[0], p[1] + epsilon * vS[1]};
		
		// choose an arbituary point a such that b=f(a) is also within 10^-6 of p
		double [] aU = divider (epsilon, p, ppU);
		double [] aS = divider (epsilon, p, ppS);
		
		double [] bU = f(aU);
		double [] bS = ffi(aS);
		
		// evaluate point b that is within 10^-6 of p
		int count = 1;
		while (distance (bS, p) > epsilon) {
			aS = divider (epsilon / ++count, p, ppS);
			bS = ffi (aS);
		}
		
		// divide line J (ab) into grids
		double [] [] aiS = new double [stable.length] [2];
		double [] [] biS = new double [stable.length] [2];
		
		aiS [0] = aS;
		biS [0] = bS;
		for (int i = 1; i < aiS.length; ++i) {
			aiS [i] = divider (i * 1.0 / stable.length, aS, bS);
			biS [i] = ffi(aiS [i]);
			
			count = 1;
			while (distance (biS [i], biS [i - 1]) > 0.001) {
				aiS[i] = divider (1.0 / ++count, aiS [i - 1], aiS [i]);
				biS[i] = ffi(aiS [i]);
			}
		}
		
		// divide the line segment ab into n grids
		double [] aiU = divider (epsilon, aU, bU);
		
		// evaluate points for the unstable manifold
		unstable [0] = f (aiU);
		for (int j = 1; j < unstable.length; ++j)
			unstable [j] = f (unstable [j - 1]);
		
		double [] pp;
		// now for each of the line J division, plot Stable manifold
		for (int i = 0; i < aiS.length; ++i) {
			stable [i] [0] = ffi (aiS [i]);
			
			for (int j = 1; j < stable [i].length; ++j) {
				pp = ffi (stable [i] [j - 1]);
				
				while (pp[0] < x0 || pp[0] > xf || pp[1] < y0 || pp[1] > yf)
					pp = ffi (pp);
				
				stable [i] [j] = pp;
			}
		}
	}
	
	// method to divide a line segment
	private double [] divider (double e, double [] p1, double [] p2)
	{
		// local variables
		double l = distance (p1, p2);
		
		double [] p = new double [p1.length];
		p[0] = p1[0] + e * (p2[0] - p1[0]) / l;
		p[1] = p1[1] + e * (p2[1] - p1[1]) / l;
		
		return p;
	}
	
	// method to calculate the distance between two points
	private double distance (double [] p1, double [] p2)
	{
		return Math.sqrt(Math.pow((p1[0] - p2[0]), 2) + Math.pow((p1[1] - p2[1]), 2));
	}
	
	// define f2(x) for negative eigenvalue functions
	private double [] ff (double [] v)
	{
		return f (f(v));
	}
	
	// define f2i(x) for negative eigenvalue function for inverse
	private double [] ffi(double[] v)
	{
		return fi(fi(v));
	}
	
	// define Henon map inverse
	private double [] fi (double [] v)
	{
		// local variable
		double [] vv = new double [2];
		vv[0] = v[1];
		vv[1] = (v[1] * v[1] + v[0] - a) / b;
		
		return vv;
	}
	
	// define Henon map
	private double [] f (double [] v)
	{
		// local variable
		double [] vv = new double [2];
		vv [0] = a - v[0] * v[0] + b * v[1];
		vv [1] = v[0];
		
		return vv;
	}
	
	// method to paint draw area on window frame
	public void paintComponent (Graphics g)
	{
		super.paintComponent(g);
		
		// show the draw area
		g.drawRect(xStart, yStart, width, height);
		g.setColor(Color.white);
		g.fillRect(xStart + 1, yStart + 1, width - 1, height - 1);
		
		// show the legend
		g.setColor(Color.black);
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
		
		// show the saddle fixed point
		g.setColor(Color.red);
		drawCross(g, (int)xP (0.94), (int)yP (0.94));
		
		// plot the unstable manifold
		g.setColor(Color.blue);
		for (int i = 0; i < unstable.length; ++i)
			g.drawOval((int)xP(unstable[i] [0]), (int) yP(unstable[i][1]), 1, 1);
		
		// plot the stable manifold
		g.setColor(Color.black);
		for (int i = 0; i < stable.length; ++i)
			for (int j = 0; j < stable [0].length; ++ j)
				g.drawOval((int)xP(stable[i] [j] [0]), (int) yP(stable[i][j][1]), 1, 1);
	}

	// method to draw a cross on the plot
	private void drawCross (Graphics g, int x, int y)
	{
		int mark = 10;
		g.drawLine(x - mark, y, x + mark, y);
		g.drawLine(x, y - mark, x, y + mark);
	}
	
	// method to convert xValue to xPixel
	private double xP (double x)
	{
		return xStart + (x - x0) / (xf - x0) * width;
	}
	
	// methdo to convert yValue to yPixel
	private double yP (double y)
	{
		return yStart + (y - yf) / (y0 - yf) * height;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		SwingUtilities.invokeLater(new Runnable () {
			public void run () {
				JFrame frame = new JFrame ("Stable and Unstable Manifold of Henon Map where " +
						"a = 2.12, b = -0.3");
				
				// instantiate an object of this
				HenonMap plot = new HenonMap (-2.5, 2.5, -2.5, 2.5);
				frame.add(plot);
				frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
				frame.setSize(620, 600);
				frame.setLocationRelativeTo(null);
				frame.setVisible(true);
			}
		});
	}

}
