package polarSystem;

import javax.swing.*;

import java.awt.*;
import java.awt.event.*;
import java.text.DecimalFormat;

public class wPolarSystem extends JPanel 
	implements MouseListener, MouseMotionListener
{
	// instance variables
	private static double y0 = -10.0;
	private static double yf = 10.0;
	private static double x0 = -10.0;
	private static double xf = 10.0;
	
	private int xStart = 50;
	private int yStart = 20;
	private int height = 500;
	private int width = 500;
	private static int screenWidth = 630;
	private static int screenHeight = 610;
	
	private double [] t;
	private double [] x;
	private double [] y;
	
	int n = 30000;
	private double dt = 0.01;
	boolean multi_plot = false;
	double a, b;
	double [] v0;
	
	private static JLabel mouseLabel = new JLabel ("");
	int fileCount = 1;
	
	// constructor
	public wPolarSystem (boolean multi, 
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
		if (v0 == null || a == 0 || b == 0)
			return;
		
		double v [] = new double [v0.length];
		v = v0;
		
		for (int i = 1; i < x.length; ++i) {
			t [i] = t [i - 1] + dt;
			v = rungeKutta (v, t [i - 1]);
			
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
		
		// the polar system
		vv[0] = v[0] * (a - v[0]);
		vv[1] = Math.pow(Math.sin(b * v[1]), 2) + 
				Math.pow((v[0] - a), 2);
		
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
			
			// draw y-axis
			g.drawLine(
					(int) (xStart + -x0 / (xf - x0) * width), yStart + 10,
					(int) (xStart + -x0 / (xf - x0) * width), 
					yStart + height - 10);
			
			// draw x-axis
			g.drawString("" + xf, xStart + width - 15, yStart + height + 15);
			g.drawLine(xStart + 10, 
					(int) (yStart + yf / (yf - y0) * height),
					xStart + width - 10, 
					(int) (yStart + yf / (yf - y0) * height));
		} else
		
		// draw plot on the panel
		for (int i = 1; i < y.length; ++i){
			// calculate the corresponding cartesian coordinates
			double xx = x[i] * Math.cos(y[i]);
			double yy = x[i] * Math.sin(y[i]);
			
			if (yy <= yf && yy >= y0 && xx >= x0 && xx <= xf)
				g.drawLine (
					(int) (xStart + (x[i - 1] * Math.cos(y[i - 1]) - x0) / (xf - x0) * width),
					(int) (yStart + (x[i - 1] * Math.sin(y[i - 1]) - yf) / (y0 - yf) * height),
					(int) (xStart + (x[i] * Math.cos(y[i]) - x0) / (xf -x0) * width),
					(int) (yStart + (x[i] * Math.sin(y[i]) - yf) / (y0 - yf) * height));
		}
	}
public void mouseDragged(MouseEvent e) {}
	
	public void mouseMoved(MouseEvent e) {
		// TODO Auto-generated method stub
		setMouseLabel (e.getX(), e.getY());
	}

	private void setMouseLabel (int x, int y)
	{
		if (x >= xStart && x <= (xStart + width) &&
				y >= yStart && y <= (yStart + height) ) {
			double xx = x0 + (xf - x0) * (x - xStart) / width;
			double yy = yf + (y0 - yf) * (y - yStart) / height;
			
			DecimalFormat f = new DecimalFormat ("#.##");
			mouseLabel.setText("(" + 
					f.format(xx) + ", " + f.format(yy) + ")");
		}
	}
	
	public void mouseClicked(MouseEvent e) {
		// TODO Auto-generated method stub
		double xp = e.getX();
		double yp = e.getY();
		
		// initialize x- and y- coordinate
		double xx = x0 + (xf - x0) * (xp - xStart) / width;
		double yy = yf + (y0 - yf) * (yp - yStart) / height;
		
		// convert to polar coordinate
		x[0] = Math.sqrt(xx * xx + yy * yy);
		if (xx > 0)
			if (yy >= 0)
				y[0] = Math.atan(yy/xx);
			else
				y[0] = 2.0 * Math.PI - Math.atan(-yy/xx);
		else if (xx < 0)
			if (yy >= 0)
				y[0] = Math.PI - Math.atan(yy/-xx);
			else
				y[0] = Math.PI + Math.atan(yy/xx);
		else
			if (yy >= 0)
				y[0] = Math.PI / 2.0;
			else
				y[0] = 3.0 * Math.PI / 2.0;
		multi_plot = true;
		v0 = new double [] {x[0], y[0]};
		solveEquation (v0);
	}

	public void mousePressed(MouseEvent e) {}

	public void mouseReleased(MouseEvent e) {}

	public void mouseEntered(MouseEvent e) {}
	
	public void mouseExited(MouseEvent e) {}
	
	// method to launch the application
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		JFrame frame = new JFrame (
				"The Phase Plot of van der Pol Oscillator");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		wPolarSystem system = new wPolarSystem (false, x0, xf, y0, yf);
		frame.add(system, BorderLayout.CENTER);
		frame.add(mouseLabel, BorderLayout.SOUTH);
		frame.setSize (screenWidth, screenHeight);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
		
		// create the controlling frame
		JFrame control = new ControlFrame (
				"Nullcline Vector Field Control", system);
		control.setSize(250, 180);
		control.setVisible (true);
	}
	
	// the control Frame
	private static class ControlFrame extends JFrame
	{
		private wPolarSystem mField;
		
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
		public ControlFrame (String title, wPolarSystem vectorField)
		{
			super (title);
			mField = vectorField;
			initComponents ();
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

	        xField.addActionListener(new ActionListener () {

				public void actionPerformed(ActionEvent e) {
					// TODO Auto-generated method stub
					processInput (e);
				}
	        	
	        });
	        yField.addActionListener(new ActionListener () {

				public void actionPerformed(ActionEvent e) {
					// TODO Auto-generated method stub
					processInput (e);
				}
	        	
	        });
	        
	        exitButton.setText("Clear");
	        exitButton.addActionListener(new java.awt.event.ActionListener() {
	            public void actionPerformed(java.awt.event.ActionEvent evt) {
	                clearButtonActionPerformed(evt);
	            }
	        });

	        drawButton.setText("Draw");
	        drawButton.addActionListener(new java.awt.event.ActionListener() {
	            public void actionPerformed(java.awt.event.ActionEvent evt) {
	                processInput(evt);
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

	    private void clearButtonActionPerformed(java.awt.event.ActionEvent evt) {                                           
	        // TODO add your handling code here:
	    	mField.multi_plot = false;
	    	mField.repaint ();
	    }                                          

	    private void processInput(java.awt.event.ActionEvent evt) {                                           
	        // TODO add your handling code here:
	    	// accept user's input
	    	try {
	    		mField.a = Double.parseDouble(xField.getText());
	    		mField.b = Double.parseDouble(yField.getText());
	    	} catch (NumberFormatException e) {
	    		JOptionPane.showMessageDialog(this, e.getMessage(),
	    				"Error Message", JOptionPane.ERROR_MESSAGE);
	    	}
	    }                                          

	}
}