package van_der_pol;

import javax.swing.*;

import java.awt.*;
import java.util.concurrent.Executors;

public class VectorField extends JPanel 
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
	
	private int [] [] x;
	private int [] [] y;

	int xGrid;
	int yGrid;
	boolean drawArrow = false;
	private double c;
	
	// constructor
	public VectorField (double x0, double xf, double y0, double yf)
	{
		// instantiate the instance variables
		this.x0 = x0;
		this.y0 = y0;
		this.xf = xf;
		this.yf = yf;
	}
	
	// method to solve the trajectory
	public void solveTrajectory ()
	{
		double [] s = new double [2];
		double d;
		
		x = new int [xGrid] [yGrid];
		y = new int [xGrid] [yGrid];
		
		c = 0.5 * Math.pow((
				Math.pow((xf - x0)/xGrid, 2) +
				Math.pow((yf - y0)/yGrid, 2)), 0.5);
		// run along the x-axis
		for (int i = 1; i < xGrid; ++i) {
			double xx = x0 + i * (xf - x0) / xGrid;
			
			// run along the y-axis
			for (int j = 1; j < yGrid; ++j) {
				double yy = y0 + j * (yf - y0) / yGrid;
				
				// evaluate the slope
				s = slope (new double [] {xx, yy});
				d = Math.pow((s[0] * s[0] + s[1] * s[1]), 0.5);
				double mod = (d ==0 ? 0 : (c / d));
				
				x[i][j] = xStart + (int)(
						width * (xx + mod * s[0] - x0) / (xf -x0));
				y[i][j] = yStart + height - (int)(
						height * (yy + mod * s[1] - y0) / (yf - y0));
			}
		}
		repaint ();
	}
	
	// equation definition
	private double [] slope (double [] v)
	{
		int l = v.length;
		double [] vv = new double [l];
		
		vv[0] = v[1];
		vv[1] = -v[0] + (1 - v[0] * v[0]) * v[1];
		
		return vv;
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
		
		// draw plot on the panel
		for (int i = 1; i < xGrid; ++i) {
			int xx = xStart + (int) (i * width / xGrid);
			
			for (int j = 1; j < yGrid; ++j) {
				int yy = yStart + height - height * j / yGrid;
				
				// draw the nullcline
				if (drawArrow)
					drawArrow(g, xx, yy, x[i][j], y[i][j]);
				else
					g.drawLine(xx, yy, x[i][j], y[i][j]);
			}
		}
	}
	
	// method to draw directional arrow
	private void drawArrow (Graphics g, int x1, int y1, int x2, int y2)
	{
		// define the array of the points that define the arrow head
		int [] xPoints = new int [3];
		int [] yPoints = new int [3];
		double theta;					// orientation of the arrow
		double gamma;					// the angle of the arrow head
		
		int xc, x3, x4, yc, y3, y4;		// the x- and y-coordinates
		double d;						// drawing division i.e sixth of arrowlength
		
		// assign values 
		d = Math.sqrt(
				Math.pow((x2 - x1), 2) + Math.pow((y2 - y1), 2)) / 6.0;
		
		if (x2 > x1) {
			if (y2 >= y1)
				 theta = Math.atan((y2 - y1) / (x2 - x1));
			else
				theta = 2.0 * Math.PI - Math.atan((y1 - y2) / (x2 - x1));
		} else if (x2 < x1) {
			if (y2 >= y1)
				theta = Math.PI - Math.atan((y2 - y1) / (x1 - x2));
			else
				theta = Math.PI + Math.atan((y1 - y2) / (x1 - x2));
		} else {
			if (y2 >= y1)
				theta = Math.PI / 2.0;
			else
				theta = 3.0 * Math.PI / 2.0;
		}
		
		gamma = Math.PI / 2.0 - theta;
		xc = (int) ((x1 + 2 * x2) / 3.0);
		yc = (int) ((y1 + 2 * y2) / 3.0);
		
		x3 = (int) (xc - 0.5 * d * Math.cos(gamma));
		y3 = (int) (yc + 0.5 * d * Math.sin(gamma));
		
		x4 = (int) (xc + 0.5 * d * Math.cos(gamma));
		y4 = (int) (yc - 0.5 * d * Math.sin(gamma));
		
		xPoints [0] = x2;
		xPoints [1] = x3;
		xPoints [2] = x4;
		
		yPoints [0] = y2;
		yPoints [1] = y3;
		yPoints [2] = y4;
		
		g.fillPolygon(xPoints, yPoints, xPoints.length);
		g.drawLine(x1, y1, x2, y2);
	}
	
	// method to launch the application
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		JFrame frame = new JFrame (
				"VectorField");
		VectorField field = new VectorField (x0, xf, y0, yf);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(field, BorderLayout.CENTER);
		
		frame.setSize (screenWidth, screenHeight);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
		
		// create the controlling frame
		JFrame control = new ControlFrame (
				"Nullcline Vector Field Control", field);
		control.setSize(250, 180);
		control.setVisible (true);
	}
	
	// the control Frame
	private static class ControlFrame extends JFrame
	{
		private VectorField mField;
		
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
		public ControlFrame (String title, VectorField vectorField)
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

	        jLabel1.setText("xGrid: ");

	        arrowCheck.setText("Show Arrow Direction");
	        arrowCheck.setHorizontalTextPosition(javax.swing.SwingConstants.LEADING);

	        jLabel3.setText("yGrid: ");

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
	    		mField.xGrid = Integer.parseInt(xField.getText());
	    		mField.yGrid = Integer.parseInt(yField.getText());
	    		mField.drawArrow = arrowCheck.isSelected();
	    		mField.solveTrajectory();
	    	} catch (NumberFormatException e) {
	    		JOptionPane.showMessageDialog(this, e.getMessage(),
	    				"Error Message", JOptionPane.ERROR_MESSAGE);
	    	}
	    }                                          

	}
}