package math.tools;

/**
 * This is the representation of a complex number with real and imaginary
 * values. This class entails static methods which are used for
 * performing operations on complex numbers
 * 
 * @author Oluwafemi Nehemiah
 * Rotneil IT Consult
 * rotneil@yahoo.com
 *
 */
public class ComplexNumber extends java.lang.Number 
implements Comparable<ComplexNumber>
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private double re, im;
	
	/**
	 * This is the no-argument constructor for initialing complex
	 * numbers with zero-real and imaginary values
	 */
	public ComplexNumber() {
		// TODO Auto-generated constructor stub
		this (0.0, 0.0);
	}
	
	/**
	 * This constructor instantiates a complex number that has only
	 * the real number and no imaginary value.
	 * 
	 * @param real The value of the real number
	 */
	public ComplexNumber (double real) {
		this (real, 0);
	}
	
	/**
	 * This is the full implementation of the ComplexNumber object.
	 * It has the real and the imaginary values
	 * @param real
	 * @param imaginary
	 */
	public ComplexNumber (double real, double imaginary) {
		re = real;
		im = imaginary;
	}
	
	// GET METHODS
	public double getRealNumber () {return this.re; }
	public double getImaginaryNumber() { return this.im; }
	
	/**
	 * This method returns the modulus of the complex number
	 * @return
	 */
	public double getModulus ()
	{
		return Math.pow((re * re + im * im), 0.5);
	}
	
	/**
	 * This method returns the conjugate of a complex matrix
	 * @param c
	 * @return
	 */
	public ComplexNumber getConjugate ()
	{
		return new ComplexNumber (getRealNumber(), -getImaginaryNumber());
	}
	
	/**
	 * This method returns the angle of the complex number as
	 * c = x + iy where theta = arctan {y/x}
	 * @return
	 */
	public double getAngleTheta ()
	{
		return Math.atan2(im, re);
	}
	
	/**
	 * This method computes angle theta in
	 * @param c
	 * @return
	 */
	public static double getAngleTheta (ComplexNumber c)
	{
		return c.getAngleTheta();
	}
	
	/**
	 * This is a utility method that returns the complex conjugate 
	 * of a complex number
	 * @param c
	 * @return
	 */
	public static ComplexNumber getConjugate(ComplexNumber c)
	{
		return c.getConjugate();
	}
	
	/**
	 * This method checks whether the complex number argument is a conjugate
	 * of the current number
	 * @param c
	 * @return
	 */
	public boolean isConjugate (ComplexNumber c)
	{
		if (c.getRealNumber() == re && c.getImaginaryNumber() == -im)
			return true;
		return false;
	}
	
	/**
	 * Checks whether c1 and c2 are conjugates of each other
	 * @param c1
	 * @param c2
	 * @return
	 */
	public static boolean areConjugate (ComplexNumber c1, ComplexNumber c2)
	{
		return c1.isConjugate(c2);
	}
	
	/**
	 * This method returns the addition of ComplexNumber this and operand,
	 * i.e this + operand.
	 * 
	 * @param operand
	 * @return
	 */
	public ComplexNumber add (ComplexNumber operand)
	{
		return new ComplexNumber (
				re + operand.getRealNumber(), 
				im + operand.getImaginaryNumber());
	}
	
	/**
	 * This method returns the addition of a double with this complex number
	 * @param operand
	 * @return
	 */
	public ComplexNumber add (double operand)
	{
		return new ComplexNumber (operand + re, im);
	}
	
	/**
	 * This methods subtracts ComplexNumber operand from this and returns the
	 * result. i.e this - operand = result
	 * @param operand
	 * @return
	 */
	public ComplexNumber subtract (ComplexNumber operand)
	{
		return new ComplexNumber (
				re - operand.getRealNumber(),
				im - operand.getImaginaryNumber()
		);
	}
	
	/**
	 * This method subtracts a double value from the real number of a complex
	 * number. result = this - value
	 * @param operand
	 * @return
	 */
	public ComplexNumber subtract (double operand)
	{
		return new ComplexNumber (re - operand, im);
	}
	
	/**
	 * This method returns the result of dividing complex number this by
	 * operand. i.e result = this / operand.
	 * 
	 * @param operand
	 * @return
	 */
	public ComplexNumber divide (ComplexNumber operand)
	{
		ComplexNumber complement = operand.getConjugate();
		ComplexNumber product = product(complement);
		double re = operand.getRealNumber(), im = operand.getImaginaryNumber();
		
		return product.divide((re * re + im * im));
	}
	
	/**
	 * This method returns the result of dividing complex number this by
	 * operand. i.e result = this / operand.
	 * 
	 * @param operand
	 * @return
	 */
	public ComplexNumber divide (double operand)
	{
		return new ComplexNumber (re / operand, im / operand);
	}
	
	/**
	 * This method returns the product of this complex number with
	 * another complex number c
	 * @param c
	 * @return
	 */
	public ComplexNumber product (ComplexNumber c)
	{
		return new ComplexNumber (
				re * c.getRealNumber() - im * c.getImaginaryNumber(),
				re * c.getImaginaryNumber() + im * c.getRealNumber()
		);
	}
	
	/**
	 * This method returns the product of a ComplexNumber with a double
	 * @param value
	 * @return
	 */
	public ComplexNumber product (double value)
	{
		return new ComplexNumber (value * re, value * im);
	}
	
	/**
	 * This method returns the product of a ComplexNumber with a double
	 * @param value
	 * @return
	 */
	public static ComplexNumber productOf (ComplexNumber c1, ComplexNumber c2)
	{
		return c1.product(c2);
	}
	
	/**
	 * This method computes the exponential of a complex number using the
	 * argument that a ComplexNumber z = x + iy can be written as Z = re(i0)
	 * so that z^n = r^n * e (in0)
	 * 
	 * @param index
	 * @return
	 */
	public ComplexNumber pow (int index)
	{
		if (index == 0.0)
			return new ComplexNumber (1.0, 0.0);
		
		ComplexNumber exp = new ComplexNumber (re, im);
		for (int i = 1; i < index; ++i)
			exp = exp.product(this);
		return exp;
	}
	
	/**
	 * This method computes the exponential of a complex number using the
	 * argument that a ComplexNumber z = x + iy can be written as Z = re(i0)
	 * so that z^n = r^n * e (in0)
	 * 
	 * @param index
	 * @return
	 */
	public ComplexNumber pow (double index)
	{
		if (index == 0.0)
			return new ComplexNumber (1.0, 0.0);
		
		double re = Math.pow(getModulus (), index);
		double theta = getAngleTheta () * index;
		
		return new ComplexNumber (re * Math.cos(theta), re * Math.sin(theta));
	}
	
	@Override
	public int intValue() {
		// TODO Auto-generated method stub
		return (int) getModulus();
	}

	@Override
	public long longValue() {
		// TODO Auto-generated method stub
		return (long) getModulus();
	}

	@Override
	public float floatValue() {
		// TODO Auto-generated method stub
		return (float) getModulus();
	}

	@Override
	public double doubleValue() {
		// TODO Auto-generated method stub
		return getModulus();
	}
	
	@Override
	public boolean equals (Object o)
	{
		ComplexNumber complexObj = (ComplexNumber) o;
		if (complexObj.getRealNumber() == re &&
				complexObj.getImaginaryNumber() == im)
			return true;
		return false;
	}
	
	public boolean equals (double value)
	{
		if (re == value && im == 0)
			return true;
		return false;
	}

	@Override
	public int compareTo(ComplexNumber o) {
		// TODO Auto-generated method stub
		// cast o to ComplexNumber
		double mod = getModulus();
		double compareMod = o.getModulus();
		if (mod < compareMod)
			return -1;
		else if (mod == compareMod)
			return 0;
		else
			return 1;
	}
	
	/**
	 * This method returns a formated ComplexNumber toString
	 */
	public String format () {
		if (im > 0)
			return String.format("%f + i%f", re, im);
		else if (im < 0) 
			return String.format("%f - i%f", re, Math.abs(im));
		else
			return String.format("%f", re);
	}
	
	@Override
	public String toString () {
		if (im > 0)
			return re + " + i" + im;
		else if (im < 0)
			return re + " - i" + Math.abs(im);
		else
			return "" + re;
	}
}
