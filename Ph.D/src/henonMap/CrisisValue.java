/*
 * This Program is used to graph the distance between the two pieces of 
 * the henon attractor with value of a around the crisis value.
 */

package henonMap;

import javax.swing.JOptionPane;

public class CrisisValue extends javax.swing.JFrame
{
	// instance variables
	private int xStart = 50;
	private int yStart = 20;
	private int width = 500;
	private int height = 500;
	
	private double x0 = -2.5;
	private double xf = 2.5;
	private double y0 = x0;
	private double yf = xf;
	
	private double a = 2.00;
	private double b = -0.3;
	private double p = 0.0;
	
	private int n = 1000000;
	private double [] [] v = new double [n][2];
	private double xx0, xxf, yy0, yyf;
	private boolean drawLine = false;
	
	public CrisisValue() {
		// set the new value of a
		setValueOfA ();
		
		// set the attractor
		setAttractor ();
		
		// create the frame
		initComponents ();
		setLocationRelativeTo (null);
	}
	
	// method to set a new value for a
	private void setValueOfA ()
	{
		try {
			this.a = Double.parseDouble(JOptionPane.showInputDialog("Enter the value of a"));
			p = (-1.3 + Math.sqrt(sq(1.3) + 4 * a)) / 2.0;
			setTitle ("The chaotic attractor of Henon Map for a = " + a + " and b = " + b);
			
		} catch (IllegalArgumentException e) {}
	}
	
	private void setAttractor ()
	{
		v [0] = new double [] {0.91, 0.91};
		for (int i = 1; i < v.length; ++i)
			v[i] = f(v[i - 1]);
	}
	
	// define the Ikeda orbit
	private double [] f (double [] v)
	{
		// local variable
		double [] vv = new double [2];
		vv [0] = a - sq(v[0]) + b * v[1];
		vv [1] = v[0];
		
		return vv;
	}
	
	// method to calculate the distance between two points
	private double distance (double [] p1, double [] p2)
	{
		double sum = 0.0;
		for (int i = 0; i < p1.length; ++i)
			sum += sq(p1[i] - p2[i]);
		
		return Math.sqrt(sum);
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
			
			// draw the saddle fixed point
			g.setColor(java.awt.Color.red);
			drawCross (g, (int)xP(p), (int)yP(p));
			
			// plot the attractor
			g.setColor(java.awt.Color.black);
			for (int i = 0; i < v.length; ++i)
				g.drawOval((int) xP(v[i][0]), (int) yP(v[i][1]), 1, 1);
			
			if (drawLine) {
				g.setColor (java.awt.Color.blue);
				g.drawLine((int) xx0, (int) yy0, (int) xxf, (int) yyf);
				distanceLabel.setText("Distance: " + 
						(distance (new double [] {xx0, yy0}, new double [] {xxf, yyf}) / 100));
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
		newButton = new javax.swing.JButton("Set new a");
		distanceLabel = new javax.swing.JLabel ("Distance: ");
		
		progressPanel = new javax.swing.JPanel ();
		progressBar = new javax.swing.JProgressBar ();
		progressBar.setMaximum((n));
		
		setDefaultCloseOperation (javax.swing.WindowConstants.EXIT_ON_CLOSE);
		setTitle ("The chaotic attractor of Henon Map for a = " + a + " and b = " + b);
		newButton.addActionListener(new java.awt.event.ActionListener () {
			public void actionPerformed(java.awt.event.ActionEvent e) {
				// TODO Auto-generated method stub
				// request for new value for a
				drawLine = false;
				distanceLabel.setText("Distance: ");
				setValueOfA();
				setAttractor ();
				plotPanel.repaint();
			}
		});
		
		addMouseListener(new java.awt.event.MouseAdapter () {

			public void mousePressed(java.awt.event.MouseEvent e) {
				// TODO Auto-generated method stub
				drawLine = true;
				xx0 = xxf = e.getPoint().getX() - 7.5;
				yy0 = yyf = e.getPoint().getY() - 30.0;
				plotPanel.repaint();
			}
			
			public void mouseReleased(java.awt.event.MouseEvent e) {
				drawLine = false;
			}
		});
		
		addMouseMotionListener(new java.awt.event.MouseMotionAdapter () {

			public void mouseDragged(java.awt.event.MouseEvent e) {
				// TODO Auto-generated method stub
				xxf = e.getPoint().getX() - 7.5;
				yyf = e.getPoint().getY() - 30.0;
				plotPanel.repaint();
			}			
		});
		
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
				.addGroup(javax.swing.GroupLayout.Alignment.TRAILING, progressPanelLayout.createSequentialGroup()
					.addGap(30, 30, 30)
					.addComponent(distanceLabel, javax.swing.GroupLayout.DEFAULT_SIZE, 220, Short.MAX_VALUE)
					.addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
					.addComponent(newButton, javax.swing.GroupLayout.DEFAULT_SIZE, 100, Short.MAX_VALUE)
					.addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
					.addComponent(progressBar))
		);
		
		progressPanelLayout.setVerticalGroup(
            progressPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(progressPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(progressPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(progressPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(distanceLabel)
                        .addComponent(newButton))
                    .addComponent(progressBar, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
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
				new CrisisValue ().setVisible(true);
			}
		});
	}
	
	// gui components
	private MyPanel plotPanel;
	private javax.swing.JButton newButton;
	private javax.swing.JLabel distanceLabel;
	private javax.swing.JPanel progressPanel;
	private javax.swing.JProgressBar progressBar;
}
