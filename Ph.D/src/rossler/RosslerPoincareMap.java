package rossler;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.JFrame;
import javax.swing.JOptionPane;

public class RosslerPoincareMap extends javax.swing.JPanel
{
	// instance variables
	private int xStart = 55;
	private int yStart = 20;
	private int height = 450;
	private int width = 450;
	
	// plot dimension
	private double x0;
	private double xf;
	private double y0;
	private double yf;
	private int n = 1500;
	private String status = "";
	
	private double [] x, y;
	private double dt = 0.01;
	
	private final int POINCARE_SECTION = 0;
	private final int RETURN_MAP = 1;
	private int plotType = POINCARE_SECTION;
			
	// popup menue
	private javax.swing.JFrame frame;
	private javax.swing.JPopupMenu popup;
	private javax.swing.JMenuItem borderMenu;
	private javax.swing.JMenuItem angleMenu;
	private javax.swing.JMenuItem poincareMenu;
	private javax.swing.JMenuItem returnMenu;
	
	// constructor
	public RosslerPoincareMap (javax.swing.JFrame frame, 
			final double theta, double x0, double xf, double y0, double yf)
	{
		// instantiate the instance data
		setBackground (java.awt.Color.white);
		
		this.x0 = x0;
		this.xf = xf;
		this.y0 = y0;
		this.yf = yf;
		
		x = new double [n];
		y = new double [n];
		
		this.frame = frame;
		popup = new javax.swing.JPopupMenu();
		borderMenu = new javax.swing.JMenuItem("Adjust Plot Border", 'b');
		angleMenu = new javax.swing.JMenuItem("Reset Section angle", 'a');
		poincareMenu = new javax.swing.JMenuItem("Poincare Section z", 'p');
		returnMenu = new javax.swing.JMenuItem("Return Map for Rn --> Rn+1", 'r');
		
		popup.add(borderMenu);
		popup.add(angleMenu);
		popup.add(poincareMenu);
		popup.add(returnMenu);
		
		addMouseListener(new java.awt.event.MouseAdapter () {
			@Override
			public void mouseReleased (java.awt.event.MouseEvent event) {
				if (event.isPopupTrigger()) {
					popup.show(event.getComponent(), event.getX(), event.getY());
				}
			}
		});
		
		// change of section angle
		angleMenu.addActionListener(new java.awt.event.ActionListener () {
			public void actionPerformed(java.awt.event.ActionEvent e) {
				// TODO Auto-generated method stub
				// prompt for the theta
				double theta = Double.parseDouble(JOptionPane.showInputDialog(null,
						"Enter the angle of the Poincare section in degrees"));
				RosslerPoincareMap.this.frame.setTitle("Poincare Section of " + 
						"Rossler flow for theta = " + theta);
				
				// convert theta to radian
				theta = theta * Math.PI / 180.0;
				startSolution(theta);
			}
		});
		
		// change of border event
		borderMenu.addActionListener(new java.awt.event.ActionListener () {
			public void actionPerformed (java.awt.event.ActionEvent e) {
				// prompt for border
				try {
					RosslerPoincareMap.this.x0 = Double.parseDouble(JOptionPane.showInputDialog(
							RosslerPoincareMap.this.frame, "Enter x0: "));
					RosslerPoincareMap.this.xf = Double.parseDouble(JOptionPane.showInputDialog(
							RosslerPoincareMap.this.frame, "Enter xf: "));
					RosslerPoincareMap.this.y0 = Double.parseDouble(JOptionPane.showInputDialog(
							RosslerPoincareMap.this.frame, "Enter y0: "));
					RosslerPoincareMap.this.yf = Double.parseDouble(JOptionPane.showInputDialog(
							RosslerPoincareMap.this.frame, "Enter yf: "));
					repaint();
				} catch (NumberFormatException ee) {
					JOptionPane.showMessageDialog(RosslerPoincareMap.this.frame,
							ee.getClass().getSimpleName(), ee.getMessage(),
							JOptionPane.ERROR_MESSAGE);
				}
			}
		});
		
		poincareMenu.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed (java.awt.event.ActionEvent e) {
				plotType = POINCARE_SECTION;
				repaint();
			}
		});
		returnMenu.addActionListener(new java.awt.event.ActionListener () {
			public void actionPerformed (java.awt.event.ActionEvent e) {
				plotType = RETURN_MAP;
				repaint ();
			}
		});
		startSolution(theta);
	}
	
	// method to start Poincare solution
	private void startSolution (final double theta)
	{
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {
			public void run () {
				showSection (theta);
			}
		});
		executor.shutdown();
	}
	
	/**
	 * This method solves for the Poincare Section through an angle theta to the x-axis
	 * 
	 * @param theta
	 */
	private void showSection (double theta)
	{
		// pick an initial condition
		double [] v = new double [] {2.0, 2.0, 2.0};
		
		for (int i = 0; i < x.length; ++i) {
			v = getNextPoincareCrossing (theta, v);
			y[i] = v[2];
			x[i] = getRadius(v);
			
			status = "Status: " + (int)(100.0 * (i + 1) / x.length) + "%";
			
			javax.swing.SwingUtilities.invokeLater(new Runnable () {
				public void run () {
					repaint();
				}
			});
		}
	}
	
	/**
	 * This method searches for the next Poincare section crossing
	 * @param theta
	 * @param v
	 * @return
	 */
	private double [] getNextPoincareCrossing (double theta, double [] v)
	{
		// iterate until there is a crossing
		while (true) {
			v = rungeKutta (v);
			
			// get the direction vector
			double [] vt = f(v);
			if (getAngleSpeed (v, vt) > 0.0 && Math.abs(theta - getAngle (v)) <= 0.0001)
				return v;
		}
	}
	
	// method to return the orbit angle in radian
	private double getAngle (double [] v)
	{
		if (v[0] == 0.0) {
			if (v[1] == 0.0)
				return 0.0;
			else
				return (v[1] > 0.0 ? Math.PI / 2.0 : 1.5 * Math.PI);
		} else if (v[0] > 0.0) {
			// first quadrant
			if (v[1] >= 0.0)
				return Math.atan(v[1] / v[0]);
			// fourth quadrant
			else
				return 2.0 * Math.PI - Math.atan(Math.abs(v[1] / v[0]));
		} else {
			// second quadrant
			if (v[1] >= 0.0)
				return Math.PI - Math.atan(Math.abs(v[1] / v[0]));
			// third quadrant
			else
				return Math.PI + Math.atan(Math.abs(v[1] / v[0]));
		}
	}
	
	// method to return the angle speed
	private double getAngleSpeed (double [] v, double [] vt)
	{
		return (v[0] * vt[1] - v[1] * vt[0]) / (v[0] * v[0] + v[1] * v[1]);
	}
	
	// method to return the radius
	private double getRadius (double [] v) 
	{
		return Math.sqrt(v[0] * v[0] + v[1] * v[1]);
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
	
	private double [] f (double [] v)
	{
		// local variables
		double a = 0.2, b = 0.2, c = 5.7;
		double [] vv = new double [v.length];
		
		vv[0] = -v[1] - v[2];
		vv[1] = v[0] + a * v[1];
		vv[2] = b + v[2] * (v[0] - c);
		
		return vv;
	}
	
	// method to paint on the panel
	public void paintComponent (java.awt.Graphics g)
	{
		super.paintComponent(g);
		
		// draw the bifurcation patterns on the screen
		if (plotType == POINCARE_SECTION) {
			g.setColor (java.awt.Color.red);
			for (int i = 0; i < x.length; ++i)
				g.drawOval(toPx(x[i]), toPy(y[i]), 1, 1);
		} else if (plotType == RETURN_MAP) {
			g.setColor (java.awt.Color.blue);
			for (int i = 1; i < x.length; ++i)
				g.drawOval(toPx(x[i - 1]), toPy(x[i]), 1, 1);
		}
		
		// show status
		g.setColor(java.awt.Color.black);
		g.drawString(status, 5, 515);
		
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
	 * This method converts x-value to the corresponding pixels
	 * @param x
	 * @return
	 */
	private int toPx (double x)
	{
		return (int) (xStart + (x - x0) / (xf - x0) * width);
	}
	
	/** 
	 * This method converts a y-value to equivalent y on the screen
	 * @param y
	 * @return
	 */
	private int toPy (double y)
	{
		return (int) (yStart + (y - yf) / (y0 - yf) * height);
	}
	
	private double toX (double Px)
	{
		return (x0 + (xf - x0) * (Px - xStart) / width);
	}
	
	private double toY (double Py)
	{
		return (yf + (y0 - yf) * (Py - yStart) / height);
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		double x0, xf, y0, yf;
		x0 = 0.0;
		xf = 10.0;
		y0 = 0.0;
		yf = 20.0;
		
		// instantiate the embedding frame
		javax.swing.JFrame frame = new javax.swing.JFrame ("Poincare Section of " + 
				"Rossler flow for theta = 0");
		
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(new RosslerPoincareMap (frame, 0.0, x0, xf, y0, yf));
		frame.setSize (580, 560);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
	}
}
