package UetaBVP;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

public class RealRegion extends JPanel {
	
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
	private int grid = 1000;
	
	private boolean [] [] xy;
	
	public RealRegion (final double x0, final double xf, 
			final double y0, final double yf)
	{
		this.x0 = x0;
		this.xf = xf;
		this.y0 = y0;
		this.yf = yf;
		
		xy = new boolean [grid] [grid];
		
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {

			public void run() 
			{
				// local variables of the grid point
				double g, k;
				
				int type = 2;
				try {
					type = Integer.parseInt(JOptionPane.showInputDialog(null,
							"Enter \n1 for 3 fixed point region\n2 for 5 fixed point region" +
							"\n3 for further period-5 fixed point\n4 for zero period-5")); 
				} catch (Exception e) {}
				
				for (int i = 0; i < grid; ++i) 
				{
					// assign the horizontal value
					g = x0 + i * (xf - x0) / grid;
					
					// assign the vertical value
					for (int j = 0; j < grid; ++j) 
					{
						k = y0 + j * (yf - y0) / grid;
						
						// set the region
						xy [i] [j] = setRegion (type, g, k);
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
	}
	
	// method to declare the basin to which an iterate has been trapped
	private boolean setRegion (int t, double g, double k)
	{
		double f5 = 120.0 / (g * k) - 95.0;
		
		if (t == 1)
			return g * k - 1 > 0;
		else if (t == 2)
			return f5 >= 0;
		else if (t == 3)
			return (5 - Math.sqrt(f5)) >= 0;
		else
			return f5 == 0.0;
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
				
				// show region
				if (xy[i][j])
					g.setColor(java.awt.Color.BLUE);
				else
					g.setColor(java.awt.Color.black);
					
				g.drawOval(xx, yy, 1, 1);
				
			}
		}
		
		// make the draw area white
		g.drawRect(xStart, yStart, width, height);
		
		// show the legend
		g.drawString("" + yf, xStart - 25, yStart + 5);
		g.drawString("" + y0, xStart - 25, yStart + height + 5);
		g.drawString("" + x0, xStart - 5, yStart + height + 15);
		g.drawString("" + xf, xStart + width - 10, yStart + height + 15);
		
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
		x0 = 0.0;
		xf = 2.0;
		y0 = 0.0;
		yf = 2.0;
		
		// TODO Auto-generated method stub
		JFrame frame = new JFrame ("The Basin of forced damped pendulum for [-pi, pi] X ["
				+ y0 + ", " + yf + "]");
		//frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(new RealRegion (x0, xf, y0, yf));
		frame.setSize (580, 560);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
	}
	
}