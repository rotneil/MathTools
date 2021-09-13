package nonlinearDyn;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Graphics;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSlider;

public class CobWebLogisticMap extends JPanel
{	
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
	
	private int period = 1;
	private double a = -1;
	
	// variables used by mouse event
	private int mx1, mx2, my1, my2;
	private boolean draw = false;
	
	public CobWebLogisticMap (final double x0, final double xf, 
			final double y0, final double yf)
	{
		this.x0 = x0;
		this.xf = xf;
		this.y0 = y0;
		this.yf = yf;
	}
	
	// define a timeT function of the system
	private double f (double x)
	{
		int count = 0;
		while (count++ < period)
			x = a - x * x;
		
		return x;
	}
	
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
		g.drawString("" + yf, xStart - 25, yStart + 5);
		g.drawString("" + y0, xStart - 25, yStart + height + 5);
		g.drawString("" + x0, xStart - 5, yStart + height + 15);
		g.drawString("" + xf, xStart + width - 10, yStart + height + 15);
		
		// show minor trends
		g.drawString("" + 1, xStart - 25, yStart + (int) (height / 4) + 5);
		g.drawString("" + -1, xStart - 25, (int) yStart + (int) (3 * height / 4) + 5);
		g.drawString("" + -1, xStart + (int) (width / 4)- 5, yStart + height + 15);
		g.drawString("" + 1, xStart + (int) (3 * width / 4)- 10, yStart + height + 15);
		
		
		// show minor marks
		for (int i = 0; i < 5; ++i) {
			// divide y-axis
			g.drawString("-", xStart, yStart + (int) (i * (height) / 4) + 4);
			g.drawString("-", xStart + width - 2, yStart + (int) (i * (height) / 4) + 4);
			g.drawString("'", xStart + (int) (i * width / 4) - 1, yStart + 10);
			g.drawString("'", xStart + (int) (i * width / 4) - 1, yStart + (height) + 7);
		}
		
		// draw the cartesian axes and the line y=x
		g.drawLine(xStart, yStart + (int) (height / 2), xStart + width, yStart + (int) (height / 2));
		g.drawLine(xStart + (int) (width / 2), yStart, xStart + (int) (width / 2), yStart + height);
		g.drawLine(xStart, yStart + height, xStart + width, yStart);
		
		// draw the function
		double x1, x2, y1, y2;
		int grid = 500;
		for (int i = 1; i < 500; ++i) {
			// assign values to x1 and x2
			x1 = x0 + (i - 1) * (xf - x0) / grid;
			x2 = x0 + i * (xf - x0) / grid;
			y1 = f (x1);
			y2 = f(x2);
			
			// ensure that the plot is within plot area
			if (y1 <= y0)
				y1 = y0;
			else if (y1 >= yf)
				y1 = yf;
			if (y2 <= y0)
				y2 = y0;
			else if (y2 >= yf)
				y2 = yf;
			g.drawLine(toXpx (x1),toYpx (y1),toXpx (x2), toYpx (y2));
		}
	}
	
	// method to convert x values to pixels
	private int toXpx (double x) 
	{
		return xStart + (int) (width * (x - x0) / (xf - x0));
	}
	
	// method to convert y values to pixels
	private int toYpx (double y) 
	{
		return yStart + (int) (height * (y - yf) / (y0 - yf));
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		double x0, xf, y0, yf;
		x0 = -2.0;
		xf = 2.0;
		y0 = -2.0;
		yf = 2.0;
		// TODO Auto-generated method stub
		JFrame frame = new JFrame ("The Cobweb plot of Logistic Map");
		
		//frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(new CobWebLogisticMap (x0, xf, y0, yf), BorderLayout.CENTER);
		
		JPanel controlPanel = new JPanel ();
		controlPanel.setLayout(new FlowLayout ());
		controlPanel.add(new JSlider (JSlider.HORIZONTAL, -1, 4, -1));
		controlPanel.add(new JLabel ("Period: "));
		frame.add(controlPanel, BorderLayout.SOUTH);
		
		frame.setSize (580, 650);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
	}
	
	// method to give a 2 s.f double
	private double todp (double x, int dp)
	{
		double sf = Math.floor(Math.pow(10, dp) * x + 0.5);
		return sf / Math.pow(10, dp);
	}
	
	private double getXValue (int xx)
	{
		return (x0 + (xf - x0) * (xx - xStart) / width);
	}
	
	private double getYValue (int yy)
	{
		return (yf + (y0 - yf) * (yy - yStart) / height);
	}
	
	
}