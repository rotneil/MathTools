package VibeResonance;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


/**
 * This program is used to compute the vibrational resonance for Duffing oscillators
 * for nonlinear double forced case
 * 		d2x(t) + d dx(t) + w0^2 x(t) + b x^3 = f cos(wt) + g cos(Wt)
 * 
 * @author Nehemiah Oluwafemi
 *
 */
public class Duffing extends javax.swing.JPanel
{
	private static final long serialVersionUID = 1L;
	// parameter declaration
	private double [] param;
	double [] Q;
	private double dt = 0.01;
	private double nT;
	
	// plot variables
	private int length = 400;
	private int start = 30;
	private double x0 = 0, xf = 500;
	private double y0 = 0, yf = 1.5;
	
	
	// constructor
	public Duffing () {
		ExecutorService executor = Executors.newCachedThreadPool();
		executor.execute(new Runnable () {
			public void run () {
				getSolution();
			}
		});
		executor.shutdown();
	}
	
	private void getSolution()
	{
		// instantiate the parameters
		param = new double [] {1.0, 1.0, 0.5, 0.1, 1.5, 15, 0};	// w0, beta, d, f, w, Omega, g
		nT = 200 * 2 * Math.PI / param[4];
		int n = (int) (nT / dt);
		
		// now generate the values of Q for values of g from zero to 500
		Q = new double [501];
		for (int i = 0; i < Q.length; ++i) {
			// evaluate the new value of g
			param [6] = i;
			double [] v = {0.1, 0.11};
			double [] x = new double [n];
			x[0] = v[0];
			
			// solve the equation
			for (int j = 1; j < x.length; ++j) {
				v = rungeKutta (v, param, dt * j);
				x[j] = v[0];
			}
			
			// now define Q
			Q[i] = qValue(x, param);
			System.out.println (i + "\t" + Q[i]);
			repaint();
		}
		
		// show the plow
		repaint ();
	}
	
	/**
	 * This method employs trapezoidal rule in order to calculate the 
	 * integer value of Q_c and Q_s
	 * 
	 * @param x The x(t)
	 * @param p The parameter values
	 * 
	 * @return
	 */
	public double qValue (double [] x, double [] p)
	{
		double vC = 0.5 * (x[0] + x[x.length - 1] * Math.cos(p[4] * dt * x.length));
		double vS = 0.5 * x[x.length - 1] * Math.sin(p[4] * dt * x.length);
		
		for (int i = 1; i < x.length - 1; ++i) {
			vC += x[i] * Math.cos(p[4] * dt * i);
			vS += x[i] * Math.sin(p[4] * dt * i);
		}
		vC *= 2.0 / nT * dt;
		vS *= 2.0 / nT * dt;
		
		return Math.sqrt(sq(vC) + sq(vS)) / param[3];
	}
	
	
	/**
     * This method uses the Runge Kutta Fourth order integeration to 
     * analyze the function defined in method f
     * 
     * @param v The start values for the scheme
     * @param p The parameters of the state function
     * @param dt
     * @return
     */
	private double [] rungeKutta (double [] v, double [] p, double t)
 	{
 		int l = v.length;	// putting in mind the time index
 		double [] c1, c2, c3, c4;
 		
 		// initialize the intermediate steps
 		c1 = new double [l];
 		c2 = new double [l];
 		c3 = new double [l];
 		c4 = new double [l];
 		
 		c1 = f(v, p, t);
 		
 		for (int i = 0; i < l; ++i)
 			c2[i] = v[i] + dt * c1[i] / 2;
 		c2 = f(c2, p, t);
 		
 		for (int i = 0; i < l; ++i)
 			c3[i] = v[i] + dt * c2[i] / 2;
 		c3 = f(c3, p, t);
 		
 		for (int i = 0; i < l; ++i)
 			c4[i] = v[i] + dt * c3[i];
 		c4 = f(c4, p, t);
 		
 		for (int i = 0; i < l; ++i)
 			c1[i] = v[i] + dt * (c1[i] + 2 * (c2[i] + c3[i]) + c4[i]) / 6.0;
 		
 		return c1;
 	}
	
	/**
	 * This is the state variable definition which may include those of the
	 * first and/or second variation equation.
	 * @param v the variables
	 * @param p The parameters
	 * @param t The time of iterate
	 * 
	 * @return
	 */
	public double [] f (double [] v, double [] p, double t)
	{
		// declare local variable
		double [] vv = new double [v.length];
		
		vv[0] = v[1];
		vv[1] = p[3] * Math.cos(p[4] * t) + p[6] * Math.cos(p[5] * t) 
				- p[2] * v[1] - sq(p[0]) * v[0] - p[1] * cube(v[0]);
		
		return vv;
	}
 	
 	private double sq (double x) {return x * x; }
	private double cube (double x) {return x * x * x;}
	
	/**
	 * This method plots the graph of Q vs g
	 */
	public void paintComponent (java.awt.Graphics g) {
		super.paintComponent(g);
		
		setBackground (java.awt.Color.white);
		
		g.drawRect(start, start, length, length);
		
		// show the legend
		g.drawString("" + yf, start - 25, start + 5);
		g.drawString("" + y0, start - 25, start + length + 5);
		g.drawString("" + x0, start - 5, start + length + 15);
		g.drawString("" + xf, start + length - 10, start + length + 15);
		
		// show the plot variables
		String xAxis = "g";
		String yAxis = "Q";
		g.drawString(xAxis, (int)(start + 0.5 * length), start + length + 20);
		g.drawString(yAxis, 15, (int)(start + 0.5 * length));
		
		// show minor marks
		for (int i = 0; i < 21; ++i) {
			// divide y-axis
			g.drawString("-", start, start + (int) (i * (length) / 20) + 4);
			g.drawString("-", start + length - 2, start + (int) (i * (length) / 20) + 4);
			g.drawString("'", start + (int) (i * length / 20) - 1, start + 10);
			g.drawString("'", start + (int) (i * length / 20) - 1, start + (length) + 7);
		}
		
		// plot the curve
		for (int i = 1; i < Q.length; ++i)
			g.drawLine(toPx(i - 1), toPy(Q[i - 1]), toPx(i), toPy(Q[i]));
	}
	
	/**
	 * This method converts x-value to the corresponding pixels
	 * @param x
	 * @return
	 */
	private int toPx (double x)
	{
		return (int) (start + (x - x0) / (xf - x0) * length);
	}
	
	/** 
	 * This method converts a y-value to equivalent y on the screen
	 * @param y
	 * @return
	 */
	private int toPy (double y)
	{
		return (int) (start + (y - yf) / (y0 - yf) * length);
	}
	/**
	 * This method launches the application
	 * @param args
	 */
	public static void main (String [] args) {
		javax.swing.JFrame frame = new javax.swing.JFrame ("Vibrational Resonance");
		frame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
		frame.add(new Duffing());
		frame.setSize(480, 500);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
	}
}
