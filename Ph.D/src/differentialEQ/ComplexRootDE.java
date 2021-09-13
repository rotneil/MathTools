package differentialEQ;

import javax.swing.*;

import java.awt.*;
import java.awt.event.*;

public class ComplexRootDE extends JPanel 
	implements MouseListener, MouseMotionListener
{
	// instance variables
	private static double y0 = -10.0;
	private static double yf = 10.0;
	private static double x0 = -3.0;
	private static double xf = 3.0;
	
	private int xStart = 50;
	private int yStart = 20;
	private int height = 450;
	private int width = 450;
	
	private double [] t;
	private double [] x;
	private double [] y;
	
	int n = 10000;
	private double dt = 0.01;
	private boolean multi_plot = false;
	
	double a, b, c, d;
	double [] v0;
	
	// constructor
	public ComplexRootDE (boolean multi, 
			double x0, double xf, double y0, double yf)
	{
		// instantiate the instance variables
		multi_plot = multi;
		this.x0 = x0;
		this.y0 = y0;
		this.xf = xf;
		this.yf = yf;
		t = new double [n];
		x = new double [n];
		y = new double [n];
		
		// register event handling
		addMouseListener (this);
		addMouseMotionListener (this);
	}
	
	// method to solve equation
	private void solveEquation (double [] v0)
	{
		double v [] = new double [v0.length];
		v = v0;
		
		for (int i = 1; i < x.length; ++i) {
			t [i] = t [i - 1] + dt;
			v = rungeKutta (v, x [i - 1]);
			
			x [i] = v [0];
			y [i] = v [1];
		}
		repaint ();
	}
	
	// Runge Kutta fourth order numerical integeration
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
	
	// equation definition
	private double [] f (double [] v, double t)
	{
		int l = v.length;
		double [] vv = new double [l];
		
		// complex roots
		vv[0] = a * v[0] - b * v[1];
		vv[1] = b * v[0] + a * v[1];
		
		return vv;
	}
	
	// method to paint on the panel
	public void paintComponent (Graphics g)
	{
		// check for multi ploting
		if (!multi_plot) {
			super.paintComponent(g);
			
			// make the draw area white
			g.setColor(Color.white);
			g.fillRect(xStart, yStart, width, height);
			g.setColor(Color.black);
			g.drawRect(xStart, yStart, width, height);
			
			// SHOW y-LEGEND
			g.drawString("" + yf, xStart - 25, (yStart + (int) (
					height * (yf - yf) / (yf - y0))));
			g.drawString("" + 0, xStart - 15, (yStart + (int) (
					height * (yf - 0) / (yf - y0))));
			g.drawString("" + y0, xStart - 25, (yStart + (int) (
					height * (yf - y0) / (yf - y0))));
			
			// SHOW x-legend
			g.drawString("" + x0, xStart - 15, yStart + height + 15);
			g.drawString("" + xf, xStart + width - 15, yStart + height + 15);
			g.drawString("" + 0, xStart + width / 2 - 5, yStart + height + 15);
			
			// show the x and y-axes
			g.drawLine(xStart + width / 2, yStart + 10,
					xStart + width / 2, yStart + height - 10);
			g.drawLine(xStart + 10, yStart + height / 2,
					xStart + width - 10, yStart + height / 2);
		}
		
		// draw plot on the panel
		for (int i = 1; i < y.length && y [i] < yf && y [i] > y0
				&& x [i] >= x0 && x [i] <= xf; ++i) {
			g.drawLine (
					(int) (xStart + (x[i - 1] - x0) / (xf - x0) * width),
					(int) (yStart + (y[i - 1] - yf) / (y0 - yf) * height),
					(int) (xStart + (x[i] - x0) / (xf -x0) * width),
					(int) (yStart + (y[i] - yf) / (y0 - yf) * height));
		}
	}
	
	public void mouseClicked(MouseEvent e) {
		// TODO Auto-generated method stub
		double xp = e.getX();
		double yp = e.getY();
		
		double xx = x0 + (xf - x0) * (xp - xStart) / width;
		double yy = yf + (y0 - yf) * (yp - yStart) / height;
		
		// initialize x- and y- coordinate
		y [0] = yy;
		x [0] = xx;
		
		multi_plot = true;
		v0 = new double [] {xx, yy};
		solveEquation (v0);
		
	}

	public void mouseDragged(MouseEvent e) {}
	public void mouseMoved(MouseEvent e) {}
	public void mousePressed(MouseEvent e) {}
	public void mouseReleased(MouseEvent e) {}
	public void mouseEntered(MouseEvent e) {}
	public void mouseExited(MouseEvent e) {}
	
	// method to launch the application
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		JFrame frame = new JFrame ("Two D phase Plot for different initial conditions");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		ComplexRootDE system = new ComplexRootDE (false, x0, xf, y0, yf);
		
		frame.add(system);
		frame.setSize (580, 560);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
		
		// create the controlling frame
		JFrame control = new ControlFrame (frame,
				"Scalar Systems of First ODE", system);
		control.setSize(250, 180);
		control.setVisible (true);
	}
	
	// the control Frame
	private static class ControlFrame extends JFrame
	{
		private ComplexRootDE mField;
		
		// gui components
		private javax.swing.JCheckBox arrowCheck;
	    private javax.swing.JButton drawButton;
	    private javax.swing.JButton exitButton;
	    private javax.swing.JLabel jLabel1;
	    private javax.swing.JLabel jLabel3;
	    private javax.swing.JPanel jPanel1;
	    private javax.swing.JPanel jPanel2;
	    private javax.swing.JTextField xField;
	    private javax.swing.JTextField yField;
		
	    // constructor
		public ControlFrame (javax.swing.JFrame frame, 
				String title, ComplexRootDE vectorField)
		{
			super (title);
			mField = vectorField;
			initComponents ();
			setLocationRelativeTo (frame);
		}
		
		private void initComponents() {

	        jPanel1 = new javax.swing.JPanel();
	        exitButton = new javax.swing.JButton();
	        drawButton = new javax.swing.JButton();
	        jPanel2 = new javax.swing.JPanel();
	        jLabel1 = new javax.swing.JLabel();
	        xField = new javax.swing.JTextField();
	        yField = new javax.swing.JTextField();
	        arrowCheck = new javax.swing.JCheckBox();
	        jLabel3 = new javax.swing.JLabel();

	        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

	        exitButton.setText("Exit");
	        exitButton.addActionListener(new java.awt.event.ActionListener() {
	            public void actionPerformed(java.awt.event.ActionEvent evt) {
	                exitButtonActionPerformed(evt);
	            }
	        });

	        drawButton.setText("Draw");
	        drawButton.addActionListener(new java.awt.event.ActionListener() {
	            public void actionPerformed(java.awt.event.ActionEvent evt) {
	                drawButtonActionPerformed(evt);
	            }
	        });

	        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
	        jPanel1.setLayout(jPanel1Layout);
	        jPanel1Layout.setHorizontalGroup(
	            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
	            .addGroup(jPanel1Layout.createSequentialGroup()
	                .addContainerGap(29, Short.MAX_VALUE)
	                .addComponent(exitButton, javax.swing.GroupLayout.PREFERRED_SIZE, 66, javax.swing.GroupLayout.PREFERRED_SIZE)
	                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
	                .addComponent(drawButton, javax.swing.GroupLayout.PREFERRED_SIZE, 67, javax.swing.GroupLayout.PREFERRED_SIZE)
	                .addGap(38, 38, 38))
	        );
	        jPanel1Layout.setVerticalGroup(
	            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
	            .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
	                .addComponent(exitButton)
	                .addComponent(drawButton))
	        );

	        jLabel1.setText("a_value: ");

	        arrowCheck.setText("Show Arrow Direction");
	        arrowCheck.setHorizontalTextPosition(javax.swing.SwingConstants.LEADING);

	        jLabel3.setText("b_value: ");

	        javax.swing.GroupLayout jPanel2Layout = new javax.swing.GroupLayout(jPanel2);
	        jPanel2.setLayout(jPanel2Layout);
	        jPanel2Layout.setHorizontalGroup(
	            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
	            .addGroup(jPanel2Layout.createSequentialGroup()
	                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
	                    .addGroup(jPanel2Layout.createSequentialGroup()
	                        .addContainerGap()
	                        .addComponent(arrowCheck, javax.swing.GroupLayout.PREFERRED_SIZE, 138, javax.swing.GroupLayout.PREFERRED_SIZE))
	                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, jPanel2Layout.createSequentialGroup()
	                        .addGap(29, 29, 29)
	                        .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
	                            .addComponent(jLabel1)
	                            .addComponent(jLabel3))
	                        .addGap(18, 18, 18)
	                        .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
	                            .addComponent(yField, javax.swing.GroupLayout.DEFAULT_SIZE, 88, Short.MAX_VALUE)
	                            .addComponent(xField))))
	                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
	        );
	        jPanel2Layout.setVerticalGroup(
	            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
	            .addGroup(jPanel2Layout.createSequentialGroup()
	                .addContainerGap()
	                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
	                    .addComponent(jLabel1)
	                    .addComponent(xField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
	                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
	                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
	                    .addComponent(yField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
	                    .addComponent(jLabel3))
	                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
	                .addComponent(arrowCheck)
	                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
	        );

	        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
	        getContentPane().setLayout(layout);
	        layout.setHorizontalGroup(
	            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
	            .addComponent(jPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
	            .addComponent(jPanel2, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
	        );
	        layout.setVerticalGroup(
	            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
	            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
	                .addComponent(jPanel2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
	                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
	                .addComponent(jPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
	        );

	        pack();
	    }// </editor-fold>                        

	    private void exitButtonActionPerformed(java.awt.event.ActionEvent evt) {                                           
	        // TODO add your handling code here:
	    	System.exit(0);
	    }                                          

	    private void drawButtonActionPerformed(java.awt.event.ActionEvent evt) {                                           
	        // TODO add your handling code here:
	    	// accept user's input
	    	try {
	    		mField.a = Integer.parseInt(xField.getText());
	    		mField.b = Integer.parseInt(yField.getText());
	    		mField.solveEquation(mField.v0);
	    	} catch (NumberFormatException e) {
	    		JOptionPane.showMessageDialog(this, e.getMessage(),
	    				"Error Message", JOptionPane.ERROR_MESSAGE);
	    	}
	    }                                          

	}
}