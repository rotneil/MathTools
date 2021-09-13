/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package UetaBVP;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.text.DecimalFormat;

/**
 *
 * @author Oluwafemi Nehemiah
 */
public class PhasePlot extends javax.swing.JPanel {
	
	// instance variables
	private int n = 100000;
	private double dt = 0.01;
	
	private double g = 1.0;
	private double k = 1.0;
	
	private double [] x;
	private double [] y;
	private java.util.ArrayList<Point> fixedPoint = new java.util.ArrayList<PhasePlot.Point>();
	
	private double x0 = -3.0;
	private double xf = 3.0;
	private double y0 = -3.0;
	private double yf = 3.0;
	
	private int xStart = 50;
	private int yStart = 50;
	private int width = 400;
	private int height = 400;
	
	private boolean multiplePlot = false;
	
    /** Creates new form UetaBVP */
    public PhasePlot () {
        // initialize the initial conditions
        x = new double [n];
        y = new double [n];
        
        // implement multiple plot here
        addMouseListener (new MouseAdapter () {
        	public void mouseClicked (MouseEvent e) {
        		// retrieve the points
        		double xx = x0 + (e.getX() - xStart) * (xf - x0) / width;
        		double yy = yf + (e.getY() - yStart) * (y0 - yf) / height;
        		
        		if (xx >= x0 && xx <= xf && yy >= y0 && yy <= yf) {
        			x[0] = xx;
        			y[0] = yy;
        			multiplePlot = true;
            		getSolution ();
        		}
        	}
        });
    }

    // method init to give the initial condition
    private void initialCondition ()
    {
    	// initial condition
    	x [0] = random (x0, xf);
    	y [0] = random (y0, yf);
    	setFixedPoint ();
    	if (fixedPoint.size() > 1)
    	{
    		x[0] = fixedPoint.get(1).getX();
    		y[0] = fixedPoint.get(1).getY();
    	} else {
    		x[0] = fixedPoint.get(0).getX();
    		y[0] = fixedPoint.get(0).getY();
    	}
    	getSolution ();
    }
    
    // method random 
    private double random (double i, double f)
    {
    	double number = Math.random();
    	number *= (f - i);
    	return (i + number);
    }
    
    // method to solve for the system and plot 
    private void getSolution ()
    {
    	double [] v = new double [] {x [0], y [0]};
    	
    	for (int i = 1; i < x.length; ++i) {
    		v = rungeKutta (v);
    		x[i] = v[0];
    		y[i] = v[1];
    	}
    	
    	// analyze the fixed points
    	setFixedPoint ();
    	analyzeFixedPoint ();
    	
    	repaint();
    }
    
    // method to analyze the fixed points
    private void analyzeFixedPoint ()
    {
    	// local variable
    	double m1, m2, m3, m4, tr, det, md;
    	String eigen;
    	analysisArea.setText("");
    	
    	for (Point p : fixedPoint) {
    		// assign value of m
    		m1 = g * Math.pow(Math.cosh(g * p.getX()), -2.0);
    		m2 = -1.0;
    		m3 = 1.0;
    		m4 = -k;
    		
    		tr = m1 + m4;
    		det = m1 * m4 - m2 * m3;
    		md = tr * tr - 4 * det;
    		
    		if (md >= 0) {
    			md = Math.sqrt(md);
    			double l1 = todp ((tr + md) / 2.0, 2);
    			eigen = "hx = " + l1 + "\t";
    			l1 = todp ((tr - md) / 2.0, 2);
    			eigen += "hy = " + l1 + "\t";
    		}
    		else {
    			md = Math.sqrt(-md);
    			double l, m, n;
    			if (tr == 0) {
    				l = todp (md / 2.0, 2);
    				eigen = "hx = i" + l + "\t";
	    			eigen += "hy = -i" + l;
    			} else {
    				l = todp (tr / 2.0, 2);
    				m = todp (md / 2.0, 2);
    				n = todp (Math.sqrt(l*l + m*m), 2);
	    			eigen = "hx = " + l + " + i" + m + "\t |d| = " + n + "\t";
	    			eigen += "hy = " + l + " - i" + m + "\t|d| = " + n;
    			}
    		}
    		analysisArea.append(p.toString() + "\t" + eigen + "\n");
    	}
    }
    
    // method to convert to decimal places
    private double todp (double value, int pow)
    {
    	int sign = value > 0 ? 1 : -1;
    	double vv = (int) (value * Math.pow(10, pow) + sign * 0.5);
    	return vv / Math.pow(10, pow);
    }
    
    // method to set fixed point
    private void setFixedPoint ()
    {
    	// the origin is a constant fixed point
    	fixedPoint.clear();
        fixedPoint.add(new Point (0.0, 0.0));
        double m = 0.0;
        
    	// the other points
        if (g * k > 1) {
        	m = Math.sqrt(3.0 * (g * k - 1.0) / (k * g * g * g));
        	
	        // use Newton-Raphson to refine the equilibrium point
	 		if (g != 0.0 && k != 0.0)
	 			m = refinePoint (g, k, m);
	 		
	 		double n = m / k;
	 		fixedPoint.add(new Point (m, n));
	 		fixedPoint.add(new Point (-m, -n));
        }
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
 	
 	// equation definition
 	private double [] f (double [] v)
 	{
 		int l = v.length;
 		double [] vv = new double [l];
 		
 		// the Ueta BVP system
 		vv[0] = -v[1] + Math.tanh(g * v[0]);
 		vv[1] = v[0] - k * v[1];
 		
 		return vv;
 	}
 	
 	// PAINT COMPONENT
 	public void paintComponent (Graphics g)
 	{
 		if (!multiplePlot) {
 			super.paintComponent(g);
 		
	 		// draw the inner rectangle
	 		g.drawRect(xStart, yStart, width, height);

	 		// draw a straight line diagonally or across the fixed points
	 		if (fixedPoint.size() > 1) {
	 			Point p1, p2;
	 			double x1, x2, y1, y2;
	 			p1 = fixedPoint.get(1);
	 			p2 = fixedPoint.get(2);
	 			x1 = p1.getX(); y1 = p1.getY();
	 			x2 = p2.getX(); y2 = p2.getY();
	 			double slope = Math.atan((y2 - y1) / (x2 - x1));
	 			int yy1 = toPy (y1 - (x1 - x0) * Math.tan(slope));
	 			int yy2 = toPy (y1 + (xf - x1) * Math.tan(slope));
	 			g.drawLine(xStart, yy1, xStart + width, yy2);
	 		}
	 		g.drawLine(xStart, yStart + height, yStart + width, yStart);
	 		
	 		// mark the equilibrium points
	 		g.setColor(Color.red);
	 		
	 		// origin
	 		for (Point p : fixedPoint)
	 			drawCross (g, p.getX(), p.getY());
 		}
 		
 		// plot the phase state
 		g.setColor(Color.black);
 		for (int i = 1; i < x.length; ++i) {
 			g.drawLine(
				(int)(xStart + (x[i - 1] - x0) / (xf - x0) * width),
				(int)(yStart + (y[i - 1] - yf) / (y0 - yf) * height),
				(int)(xStart + (x[i] - x0) / (xf - x0) * width),
				(int)(yStart + (y[i] - yf) / (y0 - yf) * height));
 		}
 	}
 	
 	// Recursive Newton-Raphson method call to 
 	private double refinePoint (double g, double k, double x0)
 	{
 		// local variable
 		if (g * k <= 1)
 			return 0.0;
 		
 		double fx = Math.tanh(g * x0) - x0 / k;
 		double ffx = g * Math.pow(Math.cosh(g * x0), -2.0) - 1.0 / k;
 		double x1 = x0 - fx / ffx;
 		
 		if (Math.abs(x1 - x0) <= 0.0001)
 			return x1;
 		else
 			return refinePoint (g, k, x1);
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
 	
 	// method to draw a cross
 	private void drawCross (Graphics g, double x, double y) 
 	{
 		// draw the horizontal line
 		g.drawLine(
			(int)(xStart + (x - 0.1 - x0) / (xf - x0) * width),
			(int)(yStart + (y - yf) / (y0 - yf) * height),
			(int)(xStart + (x + 0.1 - x0) / (xf - x0) * width),
			(int)(yStart + (y- yf) / (y0 - yf) * height));
		// draw the vertical line
		g.drawLine(
			(int)(xStart + (x - x0) / (xf - x0) * width),
			(int)(yStart + (y - 0.1 - yf) / (y0 - yf) * height),
			(int)(xStart + (x - x0) / (xf - x0) * width),
			(int)(yStart + (y + 0.1 - yf) / (y0 - yf) * height));
 		
 	}
 	
    private void kSliderStateChanged(javax.swing.event.ChangeEvent evt) {                                     
        // TODO add your handling code here:
    	// obtain the values of k- and g-
    	k = kSlider.getValue() / 100.0;
    	kLabel.setText(new DecimalFormat ("#.####").format (k));
    	multiplePlot = false;
    	initialCondition ();
    }                                    

    private void gSliderStateChanged(javax.swing.event.ChangeEvent evt) {                                     
        // TODO add your handling code here:
    	g = gSlider.getValue() / 100.0;
    	gLabel.setText(new DecimalFormat ("#.####").format(g));
    	multiplePlot = false;
    	initialCondition ();
    }                                    

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
    	try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        } catch (InstantiationException e) {
        	e.printStackTrace();
        } catch (IllegalAccessException e) {
        	e.printStackTrace();
        } catch (javax.swing.UnsupportedLookAndFeelException e) {
        	e.printStackTrace();
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
            	phasePlot = new PhasePlot ();
            	phasePlot.new MyFrame ().setVisible(true);
            }
        });
    }

    // MyFrame
    private class MyFrame extends javax.swing.JFrame {
    	public MyFrame () {
    		super ("Phase Plot of Ueta BVP System");
    		
    		jPanel2 = new javax.swing.JPanel();
            jLabel1 = new javax.swing.JLabel();
            jLabel2 = new javax.swing.JLabel();
            kSlider = new javax.swing.JSlider();
            gSlider = new javax.swing.JSlider();
            gLabel = new javax.swing.JLabel();
            kLabel = new javax.swing.JLabel();
            analysisArea = new javax.swing.JTextArea();
            jScrollPane1 = new javax.swing.JScrollPane();
            
            analysisArea.setColumns(20);
            analysisArea.setRows(5);
            jScrollPane1.setViewportView(analysisArea);

            setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
            setTitle("Phase Plot of Ueta BVP System");

            phasePlot.setBackground(new java.awt.Color(255, 255, 255));
            phasePlot.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
            phasePlot.setPreferredSize(new java.awt.Dimension(500, 500));

            javax.swing.GroupLayout phasePlotLayout = new javax.swing.GroupLayout(phasePlot);
            phasePlot.setLayout(phasePlotLayout);
            phasePlotLayout.setHorizontalGroup(
                phasePlotLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGap(0, 498, Short.MAX_VALUE)
            );
            phasePlotLayout.setVerticalGroup(
                phasePlotLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGap(0, 498, Short.MAX_VALUE)
            );

            jLabel1.setFont(new java.awt.Font("Tahoma", 1, 14)); // NOI18N
            jLabel1.setText("k");

            jLabel2.setFont(new java.awt.Font("Tahoma", 1, 14)); // NOI18N
            jLabel2.setText("g:");

            kSlider.setMajorTickSpacing(20);
            kSlider.setMaximum(200);
            kSlider.setMinorTickSpacing(2);
            kSlider.setPaintTicks(true);
            kSlider.setValue(100);
            kSlider.addChangeListener(new javax.swing.event.ChangeListener() {
                public void stateChanged(javax.swing.event.ChangeEvent evt) {
                    kSliderStateChanged(evt);
                }
            });

            gSlider.setMajorTickSpacing(20);
            gSlider.setMaximum(200);
            gSlider.setMinorTickSpacing(2);
            gSlider.setPaintTicks(true);
            gSlider.setValue(100);
            gSlider.addChangeListener(new javax.swing.event.ChangeListener() {
                public void stateChanged(javax.swing.event.ChangeEvent evt) {
                    gSliderStateChanged(evt);
                }
            });

            gLabel.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
            gLabel.setText("1.0");

            kLabel.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
            kLabel.setText("1.0");

            javax.swing.GroupLayout jPanel2Layout = new javax.swing.GroupLayout(jPanel2);
            jPanel2.setLayout(jPanel2Layout);
            jPanel2Layout.setHorizontalGroup(
                jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addComponent(jScrollPane1, javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(jPanel2Layout.createSequentialGroup()
                    .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                		.addComponent(jLabel1, javax.swing.GroupLayout.PREFERRED_SIZE, 14, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(jLabel2))
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                        .addComponent(kSlider, javax.swing.GroupLayout.DEFAULT_SIZE, 425, Short.MAX_VALUE)
                        .addComponent(gSlider, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                		.addComponent(gLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(kLabel, javax.swing.GroupLayout.DEFAULT_SIZE, 53, Short.MAX_VALUE)))
            );
            jPanel2Layout.setVerticalGroup(
                jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(jPanel2Layout.createSequentialGroup()
                    .addContainerGap()
                    .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                        .addComponent(kSlider, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(jLabel1, javax.swing.GroupLayout.PREFERRED_SIZE, 17, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(kLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 25, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGap(10, 10, 10)
                    .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                        .addComponent(gSlider, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(jLabel2)
                        .addComponent(gLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 25, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGap(0, 0, Short.MAX_VALUE)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 88, Short.MAX_VALUE)
                    )
            );

            javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
            getContentPane().setLayout(layout);
            layout.setHorizontalGroup(
                layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(layout.createSequentialGroup()
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                        .addGroup(layout.createSequentialGroup()
                            .addGap(10, 10, 10)
                            .addComponent(jPanel2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGroup(layout.createSequentialGroup()
                            .addContainerGap()
                            .addComponent(phasePlot, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                    .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            );
            layout.setVerticalGroup(
                layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(layout.createSequentialGroup()
                    .addContainerGap()
                    .addComponent(phasePlot, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                    .addComponent(jPanel2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGap(0, 0, 0))
            );

            pack();

            setLocationRelativeTo (null);
    	}	// end Frame's constructor
    	
    }	// end MyFrame
    
    // Variables declaration - do not modify                     
    private javax.swing.JSlider gSlider;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private static PhasePlot phasePlot;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JSlider kSlider;
    private javax.swing.JLabel kLabel;
    private javax.swing.JLabel gLabel;
    private javax.swing.JTextArea analysisArea;
    private javax.swing.JScrollPane jScrollPane1;
    
    // End of variables declaration 
    
    // inner class
    class Point {
    	private double x, y;
    	
    	// constructor
    	public Point (double x, double y)
    	{
    		setX (x);
    		setY (y);
    	}
    	
    	// set method
    	public void setX (double x) { this.x = x; }
    	public void setY (double y) { this.y = y; }
    	
    	// get method
    	public double getY () {return y;}
    	public double getX () {return x;}
    	
    	public String toString () {return "[" + todp(x, 2) + ", " + todp(y, 2) + "]"; }
    }

}
