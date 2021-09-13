package LorenzSystem;
import java.io.FileWriter;
import java.io.IOException;


public class Lorenz {

	public Lorenz() {
		double dt = 0.01;
		double [] v = new double [] {5.0, 5.0, 5.0};
		
		// TODO Auto-generated constructor stub
		try {
			FileWriter writer = new FileWriter (
					"C:/Users/olugbenga/Desktop/lorenz.txt");
			for (int i = 0; i < 1000; ++i) {
				double t = i * dt;
				writer.write(v [0] + "\t" + v [1] + "\t" + v [2] + "\n");
				v = rungeKutta (v, t);
			}
			writer.flush ();
			writer.close();
		} catch (IOException e) {
			System.err.println (e.getMessage());
		}
	}
	// method rungeKutta
		private double [] rungeKutta (double [] v, double t)
		{
			int l = v.length;
			double dt = 0.01;
			double [] k1 = new double [l];
			double [] k2 = new double [l];
			double [] k3 = new double [l];
			double [] k4 = new double [l];
			
			k1 = function (v, t);
			
			for (int i = 0; i < l; ++i) 
				k2 [i] = v [i] + dt * k1 [i] / 2;
			k2 = function (k2, t + dt / 2);
			
			for (int i = 0; i < l; ++i)
				k3 [i] = v [i] + dt * k2 [i] / 2;
			k3 = function (k3, t + dt / 2);
			
			for (int i = 0; i < l; ++i)
				k4 [i] = v [i] + dt * k3 [i];
			k4 = function (k4, t + dt);
			
			// the new value
			for (int i = 0; i < l; ++i)
				k1 [i] = v [i] + 
					dt * (k1 [i] + 2 * (k2 [i] + k3 [i]) + k4 [i]) / 6.0;
			
			return k1;
		}	// end method rungeKutta
		
		// define the system equation
		private double [] function (double [] v, double t)
		{
			int l = v.length;
			double [] vv = new double [l];
			vv [0] = 10 * (v [1] - v [0]);
			vv [1] = -v[0] * v [2] + 28 * v[0] - v [1];
			vv [2] = v [0] * v[1] - 8.0 * v [2] / 3.0;
			return vv;
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new Lorenz ();
	}

}
