/*
** This manifold script has removed the memory challenge of the
** previous versions of this code
*/

package Ikeda;

import javax.swing.JOptionPane;

public class Attractor extends javax.swing.JFrame
{
	// instance variables
	private int xStart = 50;
	private int yStart = 20;
	private int width = 500;
	private int height = 500;
	
	private double x0 = -1.5;//-0.75;
	private double xf = 8.0;//1.75;
	private double y0 = -4.0;//-1.5;
	private double yf = 7.5;//1.25;
	
	private int n = 1000000;
	
	private double c1 = 1.003;//0.84;
	private double c2 = 0.9;
	private double c3 = 0.4;
	private double a = 7.3;
	
	// constructor
	public Attractor ()
	{
		try {
			this.a = Double.parseDouble(JOptionPane.showInputDialog("Enter the value of a"));
		} catch (IllegalArgumentException e) {}
		
		// create the frame
		initComponents ();
		setLocationRelativeTo (null);
	}
	
	// define the Ikeda orbit
	private double [] f (double [] v)
	{
		// local variable
		double t = c3 - a / (1 + sq(v[0]) + sq(v[1]));
		double [] vv = new double [2];
		vv [0] = c1 + c2 * v[0] * Math.cos(t) - c2 * v[1] * Math.sin(t);
		vv [1] = c2 * v[0] * Math.sin(t) + c2 * v[1] * Math.cos(t);
		
		return vv;
	}
	
	private double sq (double x) { return x * x; }
	
	private class MyPanel extends javax.swing.JPanel
	{
		public void paintComponent (java.awt.Graphics g)
		{
			super.paintComponent(g);
			
			// show the draw area
			g.drawRect(xStart, yStart, width, height);
			g.setColor(java.awt.Color.white);
			g.fillRect(xStart + 1, yStart + 1, width - 1, height - 1);
			
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
			
			double v [] = new double [] {0.3, 0.4};
			for (int i = 0; i < n; ++i) {
				v = f (v);
				g.drawOval((int) xP(v[0]), (int) yP(v[1]), 1, 1);
			}
		}	// end method paintComponent
		

		// method to draw a cross on the plot
		private void drawCross (java.awt.Graphics g, int x, int y)
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
		
		// method to convert yValue to yPixel
		private double yP (double y)
		{
			return yStart + (y - yf) / (y0 - yf) * height;
		}
	}
	
	// METHOD TO CREATE THE GUI COMPONENTS
	private void initComponents ()
	{
		plotPanel = new MyPanel ();
		progressPanel = new javax.swing.JPanel ();
		progressBar = new javax.swing.JProgressBar ();
		progressBar.setMaximum((n));
		
		setDefaultCloseOperation (javax.swing.WindowConstants.EXIT_ON_CLOSE);
		setTitle ("The chaotic attractor of Ikeda Map for a = " + a);
		
		javax.swing.GroupLayout plotPanelLayout = new javax.swing.GroupLayout(plotPanel);
		plotPanel.setLayout(plotPanelLayout);
		plotPanelLayout.setHorizontalGroup(
				plotPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
				.addGap(600, 600, 600)
		);
		
		plotPanelLayout.setVerticalGroup(
				plotPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
				.addGap(540, 540, 540)
		);
		
		javax.swing.GroupLayout progressPanelLayout = new javax.swing.GroupLayout (progressPanel);
		progressPanel.setLayout(progressPanelLayout);
		progressPanelLayout.setHorizontalGroup(
				progressPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
				.addGroup (progressPanelLayout.createSequentialGroup()
					.addContainerGap(300, Short.MAX_VALUE)
					.addComponent(progressBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
					.addContainerGap())
		);
		
		progressPanelLayout.setVerticalGroup(
				progressPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
				.addGroup(progressPanelLayout.createSequentialGroup()
						.addContainerGap()
						.addComponent(progressBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
						.addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
		);
		
		javax.swing.GroupLayout layout = new javax.swing.GroupLayout (getContentPane());
		getContentPane().setLayout(layout);
		layout.setHorizontalGroup(
				layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
				.addComponent(plotPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
				.addComponent(progressPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
		);
		
		layout.setVerticalGroup(
				layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
				.addGroup(layout.createSequentialGroup()
						.addComponent(plotPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
						.addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
						.addComponent(progressPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
		);
		
		pack();
	}
	
	// METHOD TO LAUNCH APPLICATION
	public static void main (String [] args)
	{
		try {
			for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
				if (info.getName().equals("Nimbus")) {
					javax.swing.UIManager.setLookAndFeel(info.getClassName());
					break;
				}
			}
		} catch (ClassNotFoundException e) {
			
		} catch (InstantiationException e) {
			
		} catch (IllegalAccessException e) {
			
		} catch (javax.swing.UnsupportedLookAndFeelException e) {
			
		}
		
		// launch app
		java.awt.EventQueue.invokeLater(new Runnable () {
			public void run () {
				new Attractor ().setVisible(true);
			}
		});
	}
	
	// gui components
	private MyPanel plotPanel;
	private javax.swing.JPanel progressPanel;
	private javax.swing.JProgressBar progressBar;
}