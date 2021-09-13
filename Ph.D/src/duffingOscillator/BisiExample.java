package duffingOscillator;

import javax.swing.JFrame;

import sun.security.krb5.JavaxSecurityAuthKerberosAccess;

public class BisiExample extends javax.swing.JPanel {
	
	// RungeKutta solution
	public BisiExample () {
		int n = 1000000;
		double [] x = new double [n];
		double [] y = new double [n];
		double h = 0.01, w0 = 1, f = 0.2, b = 2, d = 0.1;
		double w = 1.4;
		double t = 0.0;
		
		// initial values
		x[0] = 0.1;
		y[0] = 0.11;
		
		// define the intermediate RungeKutta steps
		double k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y;
		
		// perform the n-iterations
		for (int i = 1; i < n; ++i) {
			t = t + h;
			k1x = h * y[i - 1];
			k1y = h * (f * Math.cos(w * t) - d * y[i - 1] - 
					w0 * x [i - 1] - b * Math.pow(x[i - 1], 3));
			
			k2x = h * (y[i - 1] + k1x / 2.0);
			k2y = h * (f * Math.cos(w * (t + h / 2.0)) - d * (y[i - 1] + k1y / 2.0) - 
					w0 * (x [i - 1] + k1x / 2.0) - b * Math.pow(x[i - 1] + k1x / 2.0, 3));
			
			k3x = h * (y[i - 1] + k2x / 2.0);
			k3y = h * (f * Math.cos(w * (t + h / 2.0)) - d * (y[i - 1] + k2y / 2.0) - 
					w0 * (x [i - 1] + k2x / 2.0) - b * Math.pow(x[i - 1] + k2x / 2.0, 3));
			
			k4x = h * (y[i - 1] + k3x);
			k4y = h * (f * Math.cos(w * (t + h)) - d * (y[i - 1] + k3y) - 
					w0 * (x [i - 1] + k3x) - b * Math.pow(x[i - 1] + k3x, 3));
			
			x[i] = x[i - 1] + (k1x + 2 * (k2x + k3x) + k4x) / 6.0;
			y[i] = y[i - 1] + (k1y + 2 * (k2y + k3y) + k4y) / 6.0;
		}
		/*
		for (int k = 900000; k < x.length; ++k) 
			System.out.println (x[k] + "\t" + y[k]);*/
	}
	
	public void paintComponent (java.awt.Graphics g) {
		super.paintComponent(g);
		
		setBackground (java.awt.Color.white);
		
		g.drawRect(30, 30, 400, 400);
		g.drawLine(30, 230, 430, 230);
		g.drawLine(230, 30, 230, 430);
		
		
	}
	
	public static void main (String [] args) {
		javax.swing.JFrame frame = new javax.swing.JFrame ("Duffing Oscillator");
		frame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
		frame.add(new BisiExample());
		frame.setSize(480, 500);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
		new BisiExample ();
	}
}
