package van_der_pol;

import javax.swing.*;

import java.awt.*;
import java.awt.event.*;
import java.text.DecimalFormat;

public class wLimitSet
{
	// instance variables
	private double [] t;
	private double [] x;
	private double [] y;
	
	int n = 30000;
	private double dt = 0.01;
	
	// constructor
	public wLimitSet ()
	{
		// instantiate the instance variables
		t = new double [n];
		x = new double [n];
		y = new double [n];
	}
	
	// method to solve equation
	private void solveEquation (double [] v0)
	{
		double v [] = new double [v0.length];
		v = v0;
		
		for (int i = 1; i < x.length; ++i) {
			t [i] = t [i - 1] + dt;
			v = rungeKutta (v, t [i - 1]);
			
			x [i] = v [0];
			y [i] = v [1];
		}
		
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
		
		// the Lotka-volterra equation for predatorPrey
		vv[0] = v[1];
		vv[1] = -v[0] + (1 - v[0] * v[0]) * v[1];
		
		return vv;
	}
	
	// method to launch the application
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new wLimitSet ();
	}
	
}