/*
** This manifold script has removed the memory challenge of the
** previous versions of this code
*/

package manifold;

public class ForcedDampPendulum extends javax.swing.JFrame
{
	// instance variables
	private int xStart = 50;
	private int yStart = 20;
	private int width = 500;
	private int height = 500;
	
	private double pi = Math.PI;
	
	private double x0 = -pi;
	private double xf = pi;
	private double y0 = -3.5;
	private double yf = 4.5;
	
	private double a = 2.12;
	private double b = -0.3;
	private double [] p = new double [] {-0.99, -0.33};
	private double [] vU = new double [] {1.0, -0.59};
	private double [] vS = new double [] {1.0, 0.88};
	private double epsilon = 0.000001;
	
	private int nStable = 10000;
	private int nUnstable = 10000;
	
	// constructor
	public ForcedDampPendulum ()
	{
		// create the frame
		initComponents ();
		setLocationRelativeTo (null);
	}
	
	// ethod to divide a line segment
	private double [] divider (double e, double [] p1, double [] p2)
	{
		// local variables
		double l = distance (p1, p2);
		
		double [] p = new double [p1.length];
		for (int i = 0; i < p.length; ++i)
			p[i] = p1[i] + e * (p2[i] - p1[i]) / l;
		
		return p;
	}
	
	// method to calculate the distance between two points
	private double distance (double [] p1, double [] p2)
	{
		double sum = 0.0;
		for (int i = 0; i < p1.length; ++i)
			sum += sq(p1[i] - p2[i]);
		
		return Math.sqrt(sum);
	}
	
	// method square
	private double sq (double x) {return x * x; }
	
	// define f2(x) for negative eigenvalue functions
	private double [] fft2 (double [] v)
	{
		return ft2 (ft2(v));
	}
	
	// define Henon map
	private double [] ft2 (double [] v)
	{
		// local variables
		double dt = 0.01;
		double pi2 = 2.0 * pi;
		double [] vv = new double [v.length];
		vv = v;
		
		// integrate f for the next 2-pi
		for (double t = 0.0; t <= pi2; t += dt) {
			vv = rungeKutta (vv, t, dt);
			
			// bring back theta
			int np = (int) (vv[0] / pi2 + 0.5);
			vv[0] -= pi2 * np;
		}
		
		return vv;
	}
	
	// define the inverse of f
	private double [] fit2 (double [] v)
	{
		// local variables
		double dt = -0.01;
		double pi2 = 2.0 * pi;
		double [] vv = new double [v.length];
		vv = v;
		
		// integrate f for the next 2-pi
		for (double t = 0.0; t >= -pi2; t += dt)
			vv = rungeKutta (vv, t, dt);
		
		return vv;
	}
	
	// Runge Kutta fourth order numerical integeration
	private double [] rungeKutta (double [] v, double t, double dt)
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
	
	private double [] f (double [] v, double t)
	{
		// local variable
		double [] vv = new double [2];
		/*vv [0] = v[1];
		vv [1] = 2.5 * Math.cos (t) - 0.2 * v[1] - Math.sin(v[0]);
		*/
		vv[0] = v[1];
		vv[1] = 2.5 * Math.sin(t) - 0.05 * v[1] - Math.sin(v[0]);
		
		return vv;
	}
	
	// method to paint panel background
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
			
			// show the saddle fixed point
			g.setColor(java.awt.Color.red);
			drawCross(g, (int) xP(p[0]), (int) yP(p[1]));
			
			// BEGIN THE CALCULATION OF UNSTABLE MANIFOLD
			// set the paint color to blue
			g.setColor(java.awt.Color.blue);
			
			// select a point a that is within 10^-6 of p
			double [] aU = new double [] {p[0] + epsilon * vU[0], p[1] + epsilon * vU[1]};
			aU = divider (epsilon, p, aU);
			double [] bU = ft2 (aU);
			
			// divide line ab into segment
			double [] aiU = divider (epsilon, aU, bU);
			for (int i = 0; i < nUnstable; ++i) {
				aiU = ft2 (aU);
				g.drawOval((int) xP(aiU[0]), (int) yP(aiU[1]), 1, 1);
			}
			/*
			// BEGIN THE CALCULATION OF THE STABLE MANIFOLD
			// set the paint color to black
			g.setColor(java.awt.Color.black);
			
			// select a point that is within 10^-6 of p along the stable manifold
			double [] aS = new double [] {p[0] + vS[0], p[1] + vS[1]};
			aS = divider (0.0000000001, p, aS);
			double [] bS = fit2 (aS);
			
			// divide line ab int segments
			double gridDist = distance (aS, bS) / nStable;
			double [] aiS;
			for (int i = 0; i < nStable; ++i) {
				aiS = divider (i * gridDist, aS, bS);
				
				for (int j = 0; j < 100; ++j) {
					aiS = fit2 (aiS);
					double x = aiS[0];
					double y = aiS[1];
					
					if (x >= x0 && x <= xf && y >= y0 && y <= yf)
						g.drawOval((int) xP (x), (int) yP(y), 1, 1);
				}
			}
			*/
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
		progressBar.setMaximum((nUnstable + nStable));
		
		setDefaultCloseOperation (javax.swing.WindowConstants.EXIT_ON_CLOSE);
		setTitle ("Stable and Unstable Manifold of Henon Map where a = " + a + " b = " + b);
		
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
				new ForcedDampPendulum ().setVisible(true);
			}
		});
	}
	
	// gui components
	private MyPanel plotPanel;
	private javax.swing.JPanel progressPanel;
	private javax.swing.JProgressBar progressBar;
}