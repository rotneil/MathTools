package math.tools;

/**
 * Inner abstract class for all Expression nodes
 * @author Adesola
 *
 */
public abstract class ExpNode {
	// methods that must be implemented
	abstract double value (String [] var, double [] val);
	abstract ComplexNumber complexValue (String [] var, double [] val);
	abstract StringBuffer printPostfix ();
	abstract StringBuffer printInfix ();
	abstract ExpNode derivative(String wrt);
	abstract boolean isZero ();
	abstract boolean isUnity ();
	
	String wrt () { return ""; }
	
	// utility method
	boolean hasOperator (String ex) {
		if (ex.contains("+") || ex.contains("-") || ex.contains("*")
				|| ex.contains("/") || ex.contains("%") || ex.contains("^"))
			return true;
		return false;
	}
	
	// for '+' and '-'
	boolean hasAdditionOperator (String ex) {
		if (ex.contains("+") || ex.contains("-"))
			return true;
		return false;
	}
	
	/**
	 * This method returns the parameter variables from a list of variables
	 * that contains the independent variables and the paramters.
	 * 
	 * @param x The Independent variables
	 * @param var The combination of variables
	 * 
	 * @return The parameter list.
	 */
	public static String [] getParameters (String [] x, String [] var)
	{
		// local variable
		String [] param = new String [var.length - x.length];
		int index = 0;
		for (int i = 0; i < var.length; ++i)
			if (!containsVariable (x, var[i]))
				param[index++] = var[i];
		
		return param;
	}
	
	/**
	 * A utility method used by method getParameters for checking if an 
	 * array contains a string variable str.
	 * @param array
	 * @param str
	 * @return
	 */
	public static boolean containsVariable (String [] array, String str)
	{
		for (int i = 0; i < array.length; ++i)
			if (array [i].equals(str))
				return true;
		return false;
	}
	
	/**
	 * This method returns the index of variable x in 
	 * the array of variables var.
	 * 
	 * @param x The variable whose index is sought
	 * @param var The array list of variables.
	 * @return
	 */
	public static int getIndex (String x, String [] var)
	{
		for (int i = 0; i < var.length; ++i)
			if (var[i].equals(x))
				return i;
		throw new IllegalArgumentException (
				"Variable '" + x + "' is undefined!");
	}
	
	/**
	 * This method returns the indices of the independent variables x in
	 * String array variables. This is principally used in function 
	 * integration.
	 * 
	 * @param x The independent variables
	 * @param variables The variable definition in an ExpNode
	 * @return An array of indices mapping the position of the independent
	 * variables
	 */
	public static int [] getXIndex (String [] x, String [] variables)
	{
		int [] pivot = new int [x.length];
		for (int i = 0; i < x.length; ++i)
			for (int j = 0; ; ++j)
				if (x[i].equals(variables[j])) {
					pivot [i] = j;
					break;
				}
		return pivot;
	}
	
	/**
	 * This method returns the value of a variable from variable list and values
	 * @param x The variable whose value is returned
	 * @param variables The array of variables
	 * @param values The array of a values
	 * @return The value of x
	 */
	public static double getValue (String x, String [] variables, double [] val)
	{
		for (int i = 0; i < variables.length; ++i)
			if (variables[i].equals(x))
				return val [i];
		throw new IllegalArgumentException (
				"Variable '" + x + "' is undefined!");
	}

	/**
	 * This method returns the value of an expression node represented by fx
	 * from the variable definition and values.
	 * 
	 * @param fx The ExpNode
	 * @param var Variable definition
	 * @param val Corresponding values of variables
	 * @return
	 */
	public static double getValue (ExpNode fx, String [] var, double [] val)
	{
		return fx.value(var, val);
	}
	
	/**
	 * This method returns the derivative of ExpNode fx with respect to
	 * wrt.
	 * @param fx The function to be differentiated
	 * @param wrt The independent coordinate for the differentiation
	 * @return An ExpNode containing the derivative of fx with respect to wrt
	 */
	public static ExpNode derivate (ExpNode fx, String wrt)
	{
		return fx.derivative(wrt);
	}
	
	/**
	 * The standard method that should be used in the creation of a BinOpNode.
	 * This method checks for repetition of nodes and shrinks them where 
	 * appropriate for an effective computation of values of ExpNode.
	 * 
	 * @param op The operator
	 * @param left The left ExpNode
	 * @param right The right ExpNode
	 * 
	 * @return A proper, well refined ExpNode
	 */
	public static ExpNode makeBinOpNode (char op, ExpNode left, ExpNode right)
	{
		return BinOpNode.make (op, left, right);
	}
	
	/**
	 * This method returns an ExpNode that is a variable node
	 * @param x
	 * @return
	 */
	public static VariableNode makeVariableNode (String x)
	{
		return new VariableNode (x);
	}
	
	/**
	 * This method returns an ExpNode that is a complex variable node
	 * @param re
	 * @param im
	 * @return
	 */
	public static ComplexVarNode makeComplexVarNode (String re, String im)
	{
		return new ComplexVarNode (re, im);
	}
	
	/**
	 * This method returns an ExpNode that is a constant complex number
	 * @param n
	 * @return
	 */
	public static ComplexConstantNode makeComplexConstNode (ComplexNumber n)
	{
		return new ComplexConstantNode (n);
	}
	
	/**
	 * This method returns an ExpNode that is a constant node
	 * @param a
	 * @return
	 */
	public static ConstantNode makeConstantNode (double a)
	{
		return new ConstantNode (a);
	}
	
	/**
	 * This method returns an ExpNode that is an OperandNode of the for ax^n
	 * @param a The coefficient of the Operand
	 * @param x an ExpNode
	 * @param n the index of the ExpNode
	 * @return
	 */
	public static ExpNode makeOperandNode (double a, ExpNode x, double n)
	{
		return OperandNode.make (a, x, n);
	}
	
	@Override
	public String toString () { 
		return printInfix ().toString();
	}
}


class DetNode extends ExpNode {
	// instance variable
	ExpNode [] [] X;
	
	// Constructor to form a new DetNode from an array of ExpNode
	DetNode (ExpNode [] [] x)
	{
		this.X = x;
	}
	
	/**
	 * The formation of a Determinant Node from the state equation, 
	 * independent variables x and the ComplexVar name for the 
	 * state stability eigenvalue
	 * @param fx
	 * @param x
	 * @param u
	 */
	DetNode (ExpNode [] fx, String [] x, String re, String im)
	{
		// set up the determinant X
		X = new ExpNode [fx.length] [x.length];
		ExpNode [] [] Df = Matrix.jacobian(fx, x);
		
		for (int i = 0; i < X.length; ++i)
			for (int j = 0; j < X[0].length; ++j)
				X [i] [j] = Df [i] [j];
		for (int i = 0; i < x.length; ++i)
			X [i] [i] = BinOpNode.make ('-', X[i][i], new ComplexVarNode (re, im));
	}
	
	/**
	 * This method checks if the DetNode is a zero determinant node
	 */
	@Override
	boolean isZero () {
		if (Matrix.hasZeroRow(X))
			return true;
		return false;
	}
	
	/**
	 * This method checks for ExpNode with value 1.0
	 */
	@Override
	boolean isUnity () {
		return false;
	}
	
	/**
	 * This method evaluates the determinant of this DetNode using
	 * the corresponding eigen positions to know the complex numbers.
	 * 
	 * @param var
	 * @param val
	 * @return
	 * @throws MatrixException if a problem arises in the course of 
	 * determinant computation
	 */
	@Override
	ComplexNumber complexValue (String [] var, double [] val) 
	{
		try {
			// define the complex matrix
			ComplexNumber [] [] x = new ComplexNumber [X.length] [X[0].length];
			
			// assume all the elements of the matrix are real
			for (int i = 0; i < x.length; ++i)
				for (int j = 0; j < x[0].length; ++j)
					x[i][j] = X[i][j].complexValue(var, val);
			
			// return the determinant
			return Matrix.getDeterminant(x);
			
		} catch (MatrixException e) {
			return new ComplexNumber (Double.NaN, Double.NaN);
		}
	}
	
	@Override
	double value(String[] var, double[] val) {
		// TODO Auto-generated method stub
		try {
			// get the double value of the elements
			double [] [] x = Matrix.values(X, var, val);
			return Matrix.getDeterminant(x);
		} catch (MatrixException e) {
			return Double.NaN;
		}
	}

	@Override
	StringBuffer printPostfix() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	StringBuffer printInfix() {
		// TODO Auto-generated method stub
		StringBuffer output = new StringBuffer ();
		for (int i = 0; i < X.length; ++i) {
			for (int j = 0; j < X[0].length; ++j)
				output.append(X[i][j].printInfix() + ",");
			output.append(",");
		}
		return output;
	}

	@Override
	ExpNode derivative(String wrt) {
		// TODO Auto-generated method stub
		// Express the different derivatives
		int n = X.length;
		ExpNode [] [] [] Xn = new ExpNode [n][n][X[0].length];
		
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				for (int k = 0; k < X[0].length; ++k) 
					Xn [i][j][k] = (i == j ? X[j][k].derivative(wrt) : X[j][k]);
		
		// prepare the returned result
		ExpNode result = BinOpNode.make ('+', 
				new DetNode(Xn[0]), new DetNode (Xn[1]));
		for (int i = 2; i < n; ++i)
			result = BinOpNode.make ('+', result, new DetNode(Xn[i]));
		
		return result;
	}
}

class PeriodicDetNode extends ExpNode {
	// instance variables
	ExpNode [] [] X;
	String [] x;
	String re, im;
	String param;
	
	PeriodicDetNode (String [] x, String re, String im, String param)
	{
		this.x = x;
		this.re = re;
		this.im = im;
		this.param = param;
		
		X = new ExpNode [x.length][x.length];
		
		// form the Chi function
		String [] var = Matrix.getFirstVariationIndependentVar(x);
		for (int index = x.length, i = 0; i < x.length; ++i)
			for (int j = 0; j < x.length; ++j)
				X[i][j] = new VariableNode (var[index++]);
		
		// now form the characteristic eigen value equation
		for (int i = 0; i < x.length; ++i)
			X[i][i] = BinOpNode.make ('-', X[i][i], new ComplexVarNode (re, im));
	}
	
	PeriodicDetNode (ExpNode[][] X, String[] x, String re, String im, String param) {
		this.X = X;
		this.x = x;
		this.re = re;
		this.im = im;
		this.param = param;
	}
	
	/**
	 * This method checks if the DetNode is a zero determinant node
	 */
	@Override
	boolean isZero () {
		if (Matrix.hasZeroRow(X))
			return true;
		return false;
	}
	
	/**
	 * This method checks for ExpNode with value 1.0
	 */
	@Override
	boolean isUnity () {
		return false;
	}
	
	/**
	 * This method evaluates the determinant of this DetNode using
	 * the corresponding eigen positions to know the complex numbers.
	 * 
	 * @param var
	 * @param val
	 * @return
	 * @throws MatrixException if a problem arises in the course of 
	 * determinant computation
	 */
	@Override
	ComplexNumber complexValue (String [] var, double [] val) 
	{
		try {
			// define the complex matrix
			ComplexNumber [] [] x = new ComplexNumber [X.length] [X[0].length];
			
			// initialize the complex matrix elements
			for (int i = 0; i < x.length; ++i)
				for (int j = 0; j < x[0].length; ++j)
					x[i][j] = X[i][j].complexValue(var, val);
			
			// return the determinant
			return Matrix.getDeterminant(x);
			
		} catch (MatrixException e) {
			return new ComplexNumber (Double.NaN, Double.NaN);
		}
	}
	
	@Override
	double value(String[] var, double[] val) {
		// TODO Auto-generated method stub
		try {
			// get the double value of the elements
			double [] [] x = Matrix.values(X, var, val);
			return Matrix.getDeterminant(x);
		} catch (MatrixException e) {
			return Double.NaN;
		}
	}

	@Override
	StringBuffer printPostfix() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	StringBuffer printInfix() {
		// TODO Auto-generated method stub
		StringBuffer output = new StringBuffer ();
		for (int i = 0; i < X.length; ++i) {
			for (int j = 0; j < X[0].length; ++j)
				output.append(X[i][j].printInfix() + ",");
			output.append(",");
		}
		return output;
	}

	@Override
	ExpNode derivative(String wrt) {
		// TODO Auto-generated method stub
		int n = X.length;
		ExpNode [] [] [] Xn = new ExpNode [n][n][X[0].length];
		
		// get the index of wrt in x
		// Note: method ExpNode.getIndex is not used because of the case where
		// wrt is equal to "u"
		int index = -1;
		for (int i = 0; i < x.length; ++i)
			if (x[i].equals(wrt)) {
				index = i;
				break;
			}
		
		// Express the different derivatives checking that wrt is not
		// any of the independent coordinates
		if (wrt.equals(param)) {	// differentiate wrt parameter
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					for (int k = 0; k < X[0].length; ++k)
						if (i == j) {
							if (X[j][k] instanceof VariableNode ||
									X[j][k] instanceof BinOpNode)
								Xn[i][j][k] = new VariableNode ("x" + 
								(i + 1) + "" + (k + 1) + "k");
							else
								Xn[i][j][k] = X[j][k].derivative(wrt);
						}
						else
							Xn[i][j][k] = X[j][k];
		} else if (index != -1) {	// differentiate wrt independent variable
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					for (int k = 0; k < X[0].length; ++k)
						if (i == j) {
							if (X[j][k] instanceof VariableNode ||
									X[j][k] instanceof BinOpNode)
								Xn[i][j][k] = new VariableNode ("x" + (i + 1) +
									"x" + (index < k ? (index + 1) + "" + (k + 1) :
										(k + 1) + "" + (index + 1)));
							else
								Xn[i][j][k] = X[j][k].derivative(wrt);
						}
						else
							Xn[i][j][k] = X[j][k];
		} else {	// differentiation wrt to u which can be re or im
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					for (int k = 0; k < X[0].length; ++k) 
						Xn [i][j][k] = (i == j ? X[j][k].derivative(wrt) : X[j][k]);			
		}

		// prepare the returned result
		ExpNode result = BinOpNode.make ('+', 
				new PeriodicDetNode(Xn[0], x, re, im, param), 
				new PeriodicDetNode (Xn[1], x, re, im, param));
		for (int i = 2; i < n; ++i)
			result = BinOpNode.make ('+', result, 
					new PeriodicDetNode(Xn[i], x, re, im, param));
		return result;
	}
}

class ComplexVarNode extends ExpNode {
	// an expression node that holds a number
	String re, im;
	
	// complex constructor
	ComplexVarNode (String re, String im) {
		this.re = re;
		this.im = im;
	}
	
	/**
	 * This method checks if the DetNode is a zero determinant node
	 */
	@Override
	boolean isZero () {
		return false;
	}
	
	/**
	 * This method checks for ExpNode with value 1.0
	 */
	@Override
	boolean isUnity () {
		return false;
	}
	
	@Override
	double value(String [] var, double [] val) {
		// return the modulus of the this complex number
		return complexValue (var, val).getModulus();
	}
	
	// complex number
	@Override
	ComplexNumber complexValue(String [] var, double [] val) {
		return new ComplexNumber (
				getValue (re, var, val), getValue (im, var, val));
	}
	
	@Override
	StringBuffer printPostfix() {
		return new StringBuffer (re + " + " + im);
	}

	@Override
	StringBuffer printInfix() {
		return new StringBuffer (re + " + " + im);
	}

	@Override
	ExpNode derivative(String wrt) {
		if (wrt.equals(re))
			return new ConstantNode (1.0);
		else if (wrt.equals(im))
			return new ComplexConstantNode(new ComplexNumber (0.0, 1.0));
		else
			return new ConstantNode(0.0);
	}
}

class ComplexConstantNode extends ExpNode {
	// an expression node that holds a number
	ComplexNumber number;
	
	// constructor
	ComplexConstantNode (ComplexNumber c) {
		number = c;
	}
	
	/**
	 * This method checks if the DetNode is a zero determinant node
	 */
	@Override
	boolean isZero () {
		if (number.getRealNumber() == 0.0 && number.getImaginaryNumber() == 0.0)
			return true;
		return false;
	}
	
	/**
	 * This method checks for ExpNode with value 1.0
	 */
	@Override
	boolean isUnity () {
		if (number.getRealNumber() == 1.0 && number.getImaginaryNumber() == 0.0)
			return true;
		return false;
	}
	
	@Override
	double value(String [] var, double [] val) {return number.getModulus(); }
	
	@Override
	ComplexNumber complexValue(String [] var, double [] val) { return number; }
	
	@Override
	StringBuffer printPostfix() {
		return new StringBuffer (String.valueOf(number));
	}

	@Override
	StringBuffer printInfix() {
		if (number.getRealNumber() == 0.0 
				&& number.getImaginaryNumber() == 0.0)
			return new StringBuffer("0.0");
		return new StringBuffer (number.toString());
	}

	@Override
	ExpNode derivative(String wrt) {return new ConstantNode(0.0);}
}

class ConstantNode extends ExpNode {
	// an expression node that holds a number
	double number;
	
	// constructor
	ConstantNode (double val) {
		number = val;
	}
	
	@Override boolean isZero () { return number == 0.0; }
	@Override boolean isUnity () { return number == 1.0; }
	
	@Override
	double value(String [] var, double [] val) {return number; }
	double value() { return number; }
	
	@Override
	ComplexNumber complexValue (String [] var, double [] val) {
		return new ComplexNumber (number);
	}
	
	@Override
	StringBuffer printPostfix() {
		return new StringBuffer (String.valueOf(number));
	}

	@Override
	StringBuffer printInfix() {
		if (number == 0.0)
			return new StringBuffer("0.0");
		return new StringBuffer (String.valueOf(number));
	}

	@Override
	ExpNode derivative(String wrt) {return new ConstantNode(0.0);}
}

class VariableNode extends ExpNode {
	// instance variable
	String x;
	
	// an expression node that represents a reference to the variable x
	VariableNode (String x) {
		this.x = x;
	}
	
	@Override boolean isZero () { return false; }
	@Override boolean isUnity () { return false; }
	@Override String wrt () { return x; }
	
	@Override
	double value (String [] var, double [] val) 
	{
		// get the value of the variable
		return getValue (x, var, val);
	}
	
	@Override
	ComplexNumber complexValue (String [] var, double [] val) {
		return new ComplexNumber (value(var, val));
	}
	
	@Override
	StringBuffer printPostfix() {
		return new StringBuffer (x);
	}

	@Override
	StringBuffer printInfix() {
		return new StringBuffer (x);
	}

	@Override
	ExpNode derivative(String wrt) {
		// TODO Auto-generated method stub
		if (x.equals(wrt))
			return new ConstantNode (1);
		return new ConstantNode (0);
	}
}

class OperandNode extends ExpNode {
	// instance variable
	double coef;
	double index;
	ExpNode x;
	
	// constructor that represents type ax^n
	OperandNode (double a, ExpNode x, double n) {
		this.coef = a;
		this.index = n;
		this.x = x;
	}
	
	@Override boolean isZero () { 
		if (coef == 0.0)
			return true;
		return false; 
	}
	
	@Override boolean isUnity () { return false; }
	
	/**
	 * This is a standard method which must be used to create a refined 
	 * Operand node in order to prevent polification of nodes
	 *  
	 * @param a
	 * @param x
	 * @param n
	 * @return
	 */
	static ExpNode make (double a, ExpNode x, double n) {
		// isolate the unnecessary cases
		if (a == 0.0)
			return new ConstantNode (a);
		else if (n == 0.0)
			return new ConstantNode (a);
		else if (a == 1.0 && n == 1.0)
			return x;
		else if (x instanceof ConstantNode) {
			// cast x to Constant node
			ConstantNode xNode = (ConstantNode) x;
			
			if (xNode.number == 0.0)
				return new ConstantNode (0.0);
			else if (xNode.number == 1.0)
				return new ConstantNode (a);
			else
				return new ConstantNode (a * Math.pow(xNode.number, n));
		}
		
		return new OperandNode (a, x, n);
	}
	
	@Override
	double value (String [] var, double [] val) 
	{
		// get the value of the variable
		return (coef * Math.pow(x.value(var, val), index)); 
	}
	
	@Override
	ComplexNumber complexValue (String [] var, double [] val) 
	{
		ComplexNumber c = x.complexValue(var, val);
		c = c.product(coef);
		return c.pow(index);
	}
	
	@Override
	StringBuffer printPostfix() {
		return new StringBuffer (
				(coef != 1.0 ? coef + " " : "") + x.printPostfix() + 
				(index != 1.0 ? " " + index + " ^": "") +
				(coef != 1.0 ? " *" : ""));
	}

	@Override
	StringBuffer printInfix() {
		return new StringBuffer (
				(coef != 1.0 ? coef + " * " : "") + 
				(hasOperator (x.printInfix().toString()) ? 
						"(" + x.printInfix() + ") " : x.printInfix()) +
				(index != 1.0 ? " ^ " + index : ""));
	}

	@Override
	ExpNode derivative(String wrt) {
		// TODO Auto-generated method stub
		if (index == 0)
			return new ConstantNode (0);
		if (index == 1)
			return BinOpNode.make ('*', new ConstantNode (coef), 
					x.derivative(wrt));
		return BinOpNode.make ('*', x.derivative(wrt), 
				OperandNode.make (coef * index, x, index - 1));
	}
}


class BinOpNode extends ExpNode {
	// An expression node representing a binary operator
	char op;		// The operator.
	ExpNode left;	// the expression for its left operand
	ExpNode right;	// The expression for its right operand.
	
	/**
	 * This constructor creates a Binary Operation Node for algebraic expressions.
	 * Since it does not check the left and right ExpNodes for repetition of terms,
	 * it's use is not encouraged. Instead, use the static method 
	 * @see {make (op, left, right)}
	 * 
	 * @param op The operator
	 * @param left
	 * @param right
	 */
	BinOpNode (char op, ExpNode left, ExpNode right) {
		// construct a BinOpNode containing the specified data
		this.op = op;
		this.left = left;
		this.right = right;
	}
	
	@Override
	boolean isZero () {
		switch (op) {
		case '+': case '-':
			if (left.isZero() && right.isZero())
				return true;
			break;
		case '*':
			if (left.isZero() || right.isZero())
				return true;
			break;
		case '/':
			if (left.isZero())
				return true;
			break;
		}
		return false;
	}
	
	@Override boolean isUnity () { 
		switch (op) {
		case '*':
			if (left.isUnity() && right.isUnity())
				return true;
			break;
		case '/':
			if (left.isUnity() && right.isUnity())
				return true;
			break;
		}
		return false;
	}
	
	/**
	 * The standard method that should be used in the creation of a BinOpNode.
	 * This method checks for repetition of nodes and shrinks them where 
	 * appropriate for an effective computation of values of ExpNode.
	 * 
	 * The rearrangement is done only to the first sub nodes of BinOpNodes
	 * 
	 * @param op The operator
	 * @param left The left ExpNode
	 * @param right The right ExpNode
	 * 
	 * @return A proper, refined BinOpNode
	 */
	static ExpNode make (char op, ExpNode left, ExpNode right)
	{
		switch (op) {
			case '+':
				if (left.isZero())
					return right;
				if (right.isZero())
					return left;
				break;
			case '-':
				if (left.isZero())
					return ExpNode.makeOperandNode(-1.0, right, 1.0);
				if (right.isZero())
					return left;
				break;
			case '*':
				if (left.isZero())
					return left;
				if (right.isZero() )
					return right;
				if (left.isUnity())
					return right;
				if (right.isUnity())
					return left;
				break;
			case '/':
				if (left.isZero())
					return left;
				if (right.isUnity())
					return left;
				break;
			case '^':
				if (right.isZero())
					return new ConstantNode (1.0);
				if (left.isUnity())
					return left;
				if (left.isZero())
					return left;
				break;
		}
		/*
		// treatment based on nodes of the BinOpNode
		if (left instanceof ConstantNode) {
			ConstantNode lNode = (ConstantNode) left;
			if (lNode.number == 0.0)
				switch (op) {
					case '*': case '/': return left;
					case '+': return right;
					case '-': return OperandNode.make (-1.0, right, 1.0);
				}
			
			// the case of left being 1
			if (lNode.number == 1.0 && op == '*')
				return right;
			
			if (right instanceof ConstantNode) {
				ConstantNode rNode = (ConstantNode) right;
				if (rNode.number == 0.0)
					switch (op) {
						case '*': return right;
						case '/': throw new NumberFormatException (
								"Division by zero exception");
						case '+': case '-': return left;
					}
				else if (rNode.number == 1.0 && (op == '*' || op == '/'))
					return left;
				
				// return an appropriate constant operation
				switch (op) {
					case '+': lNode.number += rNode.number; break;
					case '-': lNode.number -= rNode.number; break;
					case '*': lNode.number *= rNode.number; break;
					case '/': lNode.number /= rNode.number; break;
				}
				
				return lNode;
			} else if (right instanceof VariableNode) {
				if (op == '*')
					return OperandNode.make (lNode.number, right, 1.0);
				else if (op == '/')
					return OperandNode.make (lNode.number, right, -1.0);
			} else if (right instanceof OperandNode) {
				if (op == '*') {
					OperandNode rOp = (OperandNode) right;
					rOp.coef *= lNode.number;
					return right;
				}
			} else if (right instanceof BinOpNode) {
				if (op == '*') {
					// down cast right
					BinOpNode rightBin = (BinOpNode) right;
					if (rightBin.op == '*') {
						// work on the left node
						if (rightBin.left instanceof ConstantNode) {
							// down cast rightBin.left
							ConstantNode rightL = (ConstantNode) rightBin.left;
							
							return OperandNode.make(lNode.number * rightL.number,
									rightBin.right, 1.0);
						} else if (rightBin.right instanceof ConstantNode) {
							// down cast rightBin.left
							ConstantNode rightR = (ConstantNode) rightBin.right;
							
							return OperandNode.make(lNode.number * rightR.number,
									rightBin.left, 1.0);
						} else if (rightBin.left instanceof VariableNode) {
							rightBin.left = OperandNode.make(lNode.number, 
									rightBin.left, 1.0);
							return right;
						} else if (rightBin.right instanceof VariableNode) {
							rightBin.right = OperandNode.make(lNode.number, 
									rightBin.right, 1.0);
							return right;
						} else if (rightBin.left instanceof OperandNode) {
							// downcast
							OperandNode rightL = (OperandNode) rightBin.left;
							rightL.coef *= lNode.number;
							return right;
						} else if (rightBin.right instanceof OperandNode) {
							// downcast
							OperandNode rightR = (OperandNode) rightBin.right;
							rightR.coef *= lNode.number;
							return right;
						}
					}
				}
			}
		}
		
		// this case is the next special case in the refinement of BinOpNode.
		else if (left instanceof VariableNode) {
			// cast left node
			VariableNode lVar = (VariableNode) left;
			
			// the case in which the right is made of ConstantNode is simply
			// interchanged with left in order to reflect the lower precedence
			// refinement
			if (right instanceof ConstantNode) {
				// cast right
				ConstantNode rNode = (ConstantNode) right;
				if (rNode.number == 0.0)
					switch (op) {
					case '-': case '+': return left;
					case '*': return right;
					case '/': throw new NumberFormatException (
							"Division by zero error");
					}
				
				// treat the multiplication and division cases
				if (op == '*')
					return OperandNode.make (rNode.number, left, 1.0);
				else if (op == '/')
					return OperandNode.make ((1.0 / rNode.number), left, 1.0);
			} else if (right instanceof VariableNode) {
				// cast right node
				VariableNode rVar = (VariableNode) right;
				if (lVar.x.equals(rVar.x)) {
					switch (op) {
						case '-': return new ConstantNode (0);
						case '+': return OperandNode.make (2.0, left, 1.0);
						case '*': return OperandNode.make (1.0, left, 2.0);
						case '/': return new ConstantNode (1.0);
					}
				} 
			} else if (right instanceof OperandNode) {
				// case right node
				OperandNode rOp = (OperandNode) right;
				
				// check the ExpNode x of rOp
				if (rOp.x instanceof VariableNode) {
					// cast its ExpNode x
					VariableNode rOpX = (VariableNode) rOp.x;
					
					// check for similar Variable with left var
					if (lVar.x.equals(rOpX.x)) {
						if (op == '*')
							rOp.index += 1.0;
						else if (op == '/') {
							rOp.coef = 1.0 / rOp.coef;
							rOp.index = 1.0 - rOp.index;
						} else if (op == '+' && rOp.index == 1.0) {
							rOp.coef += 1.0;
						} else if (op == '-' && rOp.index == 1.0) {
							rOp.coef = 1.0 - rOp.coef;
						}
						return rOp;
					}
				}
			} else if (right instanceof BinOpNode) {
				// down cast right
				BinOpNode rightBin = (BinOpNode) right;
				
				if (op == '*' && rightBin.op == '*') {
					// check if any of the subleft or subright nodes is the same
					// as the left variable node
					if (rightBin.left instanceof VariableNode) {
						// down cast
						VariableNode rightLBin = (VariableNode) rightBin.left;
						if (lVar.x.equals(rightLBin.x)) {
							// set the subleft of the right as an Operand Node
							rightBin.left = new OperandNode (1.0, left, 2.0);
							return right;
						}
					} else if (rightBin.right instanceof VariableNode) {
						// down cast
						VariableNode rightRBin = (VariableNode) rightBin.right;
						if (lVar.x.equals(rightRBin.x)) {
							// set the subright of the right as an operand node
							rightBin.right = new OperandNode (1.0, left, 2.0);
							return right;
						}
					} else if (rightBin.left instanceof OperandNode) {
						// down cast
						OperandNode rightLBin = (OperandNode) rightBin.left;
						if (lVar.x.equals(rightLBin.x)) {
							// set the subleft of the right as an advanced node
							rightLBin.index += 1.0;
							return right;
						}
					} else if (rightBin.right instanceof OperandNode) {
						// down cast
						OperandNode rightRBin = (OperandNode) rightBin.right;
						if (lVar.x.equals(rightRBin.x)) {
							// set the subright of the right as an advanced node
							rightRBin.index += 1.0;
							return right;
						}
					}
				}
			}
		} 
		
		else if (left instanceof OperandNode) {
			// cast left node
			OperandNode lOpNode = (OperandNode) left;
			
			// case where the right is a ConstantNode
			if (right instanceof ConstantNode) {
				ConstantNode rNode = (ConstantNode) right;
				if (rNode.number == 0.0)
					switch (op) {
					case '-': case '+': return left;
					case '*': return right;
					case '/': throw new NumberFormatException (
							"Division by zero error.");
					}
				
				// case of a proper constant right node
				switch (op) {
				case '*':
					lOpNode.coef *= rNode.number;
					return left;
				case '/':
					lOpNode.coef /= rNode.number;
					return left;
				}
			} else if (right instanceof VariableNode) {
				// check if lOpNode is composed principally of Variable
				if (lOpNode.x instanceof VariableNode) {
					// cast lOpNode and right VariableNodes
					VariableNode leftX = (VariableNode) lOpNode.x;
					VariableNode rightX = (VariableNode) right;
					
					if (leftX.x.equals(rightX.x))
						switch (op) {
						case '*':
							lOpNode.index += 1.0;
							return lOpNode;
						case '/':
							lOpNode.index -= 1.0;
							return lOpNode;
						case '+':
							if (lOpNode.index == 1.0) {
								lOpNode.coef += 1.0;
								return lOpNode;
							}
						case '-':
							if (lOpNode.index == 1.0) {
								lOpNode.coef -= 1.0;
								return lOpNode;
							}
						}
				}
			} else if (right instanceof OperandNode) {
				// cast right
				OperandNode rOpNode = (OperandNode) right;
				if (lOpNode.x instanceof VariableNode && 
						rOpNode.x instanceof VariableNode) {
					// down cast their x-variable
					VariableNode leftVar = (VariableNode) lOpNode.x;
					VariableNode rightVar = (VariableNode) rOpNode.x;
					if (leftVar.x.equals(rightVar.x)) {
						switch (op) {
						case '*':
							lOpNode.coef *= rOpNode.coef;
							lOpNode.index += rOpNode.index;
							return left;
						case '/':
							lOpNode.coef /= rOpNode.coef;
							lOpNode.index -= rOpNode.index;
							return left;
						}
					}
				}
			}
			
		} else if (left instanceof DetNode) {
			// down cast left DetNode
			DetNode leftDet = (DetNode) left;
			
			// get the ExpNode X of the leftDet
			ExpNode [][] X = leftDet.X;
			
			if (Matrix.hasZeroRow(X)) {
				switch (op) {
					case '+' : case '-':
						return right;
					case '*':
						return new ConstantNode (0.0);
				}
			}
			
		} else if (left instanceof PeriodicDetNode) {
			// down cast left PeriodicDetNode
			PeriodicDetNode leftDet = (PeriodicDetNode) left;
			
			// get the ExpNode X of the leftDet
			ExpNode [][] X = leftDet.X;
			
			if (Matrix.hasZeroRow(X)) {
				switch (op) {
					case '+' : case '-':
						return right;
					case '*':
						return new ConstantNode (0.0);
				}
			}
			
		} else if (left instanceof BinOpNode) {
			// down cast left BinOpNode
			BinOpNode leftBin = (BinOpNode) left;
			
			// treat the right special cases
			if (right instanceof ConstantNode) {
				// downcast right
				ConstantNode rightCons = (ConstantNode) right;
				// treat the special cases with zero-right
				if (rightCons.number == 0.0)
					switch (op) {
					case '-': case '+':
						return left;
					case '*': return right;
					case '/': throw new NumberFormatException (
							"Division by zero exception");
					}
				// cases with 1.0
				if (rightCons.number == 1.0 && (op == '*' || op == '/'))
					return left;
				
				// the general constant cases
				if (op == '*') {
					if (leftBin.left instanceof ConstantNode) {
						// down cast leftBin.left
						ConstantNode leftL = (ConstantNode) leftBin.left;
						leftL.number *= rightCons.number;
						return left;
					} else if (leftBin.right instanceof ConstantNode) {
						// down cast leftBin.left
						ConstantNode leftR = (ConstantNode) leftBin.right;
						leftR.number *= rightCons.number;
						return left;
					} else if (leftBin.left instanceof VariableNode) {
						leftBin.left = OperandNode.make(rightCons.number,
								leftBin.left, 1.0);
						return left;
					} else if (leftBin.right instanceof VariableNode) {
						leftBin.right = OperandNode.make(rightCons.number,
								leftBin.right, 1.0);
						return left;
					} else if (leftBin.left instanceof OperandNode) {
						// down cast
						OperandNode leftL = (OperandNode) leftBin.left;
						leftL.coef *= rightCons.number;
						return left;
					} else if (leftBin.right instanceof OperandNode) {
						// down cast
						OperandNode leftR = (OperandNode) leftBin.right;
						leftR.coef *= rightCons.number;
						return left;
					}
				}
			} else if (right instanceof VariableNode) {
				if (op == '*' || op == '/') {
					// down cast right
					VariableNode rightVar = (VariableNode) right;
					if (leftBin.op == '*') {
						if (leftBin.left instanceof VariableNode) {
							// down cast leftBin.left
							VariableNode leftL = (VariableNode) leftBin.left;
							if (leftL.x.equals(rightVar.x)) {
								switch (op) {
								case '*':
									leftBin.left = OperandNode.make(1.0, right, 2.0);
									return left;
								case '/':
									return leftBin.right;
								}
							}
						} else if (leftBin.right instanceof VariableNode) {
							// down cast leftBin.left
							VariableNode leftR = (VariableNode) leftBin.right;
							if (leftR.x.equals(rightVar.x)) {
								switch (op) {
								case '*':
									leftBin.right = OperandNode.make(1.0, right, 2.0);
									return left;
								case '/':
									return leftBin.left;
								}
							}
						} else if (leftBin.left instanceof OperandNode) {
							// down cast leftBin.left
							OperandNode leftL = (OperandNode) leftBin.left;
							if (leftL.x instanceof VariableNode) {
								// down cast
								VariableNode leftLX = (VariableNode) leftL.x;
								if (leftLX.x.equals(rightVar.x)) {
									switch (op) {
									case '*':
										leftL.index += 1.0;
										return left;
									case '/':
										leftL.index -= 1.0;
										return left;
									}
								}
							}
						} else if (leftBin.right instanceof OperandNode) {
							// down cast leftBin.left
							OperandNode leftR = (OperandNode) leftBin.right;
							if (leftR.x instanceof VariableNode) {
								// down cast
								VariableNode leftRX = (VariableNode) leftR.x;
								if (leftRX.x.equals(rightVar.x)) {
									switch (op) {
									case '*':
										leftR.index += 1.0;
										return left;
									case '/':
										leftR.index -= 1.0;
										return left;
									}
								}
							}
						}
					}
				}
			} else if (right instanceof OperandNode) {
				// treat only * and / cases
				if (op == '*' || op == '/') {
					// downcast right
					OperandNode rightOp = (OperandNode) right;
					// check if the ExpNode in rightOp is a variable node
					if (rightOp.x instanceof VariableNode) {
						// downcast rightOp variable node
						VariableNode rightVar = (VariableNode) rightOp.x;
						
						// treat the cases of left operator being *
						if (leftBin.op == '*') {
							if (leftBin.left instanceof ConstantNode) {
								// down cast
								ConstantNode leftL = (ConstantNode) leftBin.left;
								rightOp.coef *= leftL.number;
								leftBin.left = rightOp;
								return left;
							} else if (leftBin.right instanceof ConstantNode) {
								// down cast
								ConstantNode leftR = (ConstantNode) leftBin.right;
								rightOp.coef *= leftR.number;
								leftBin.right = rightOp;
								return left;
							} else if (leftBin.left instanceof VariableNode) {
								// down cast
								VariableNode leftL = (VariableNode) leftBin.left;
								if (leftL.x.equals(rightVar.x)) {
									rightOp.index += 1.0;
									leftBin.left = rightOp;
									return left;
								}
							} else if (leftBin.right instanceof VariableNode) {
								// down cast
								VariableNode leftR = (VariableNode) leftBin.right;
								if (leftR.x.equals(rightVar.x)) {
									rightOp.index += 1.0;
									leftBin.right = rightOp;
									return left;
								}
							} else if (leftBin.left instanceof OperandNode) {
								// down cast
								OperandNode leftL = (OperandNode) leftBin.left;
								if (leftL.x instanceof VariableNode) {
									// cast it
									VariableNode leftLX = (VariableNode) leftL.x;
									if (leftLX.x.equals(rightVar.x)) {
										leftL.coef *= rightOp.coef;
										leftL.index += rightOp.index;
										return left;
									}
								}
							} else if (leftBin.right instanceof OperandNode) {
								// down cast
								OperandNode leftR = (OperandNode) leftBin.right;
								if (leftR.x instanceof VariableNode) {
									// cast it
									VariableNode leftRX = (VariableNode) leftR.x;
									if (leftRX.x.equals(rightVar.x)) {
										leftR.coef *= rightOp.coef;
										leftR.index += rightOp.index;
										return left;
									}
								}
							}
						}
					}
				}
			}
		}*/
			
		return new BinOpNode (op, left, right);
	}
	
	// allow Complex value computation
	@Override
	ComplexNumber complexValue (String [] var, double [] val) 
	{
		ComplexNumber leftComplex = left.complexValue(var, val);
		ComplexNumber rightComplex = right.complexValue(var, val);
		
		switch (op) {
		case '-':
			return leftComplex.subtract(rightComplex);
		case '+':
			return leftComplex.add(rightComplex);
		case '*':
			return leftComplex.product(rightComplex);
		case '/':
			return leftComplex.divide(rightComplex);
		default:
			return new ComplexNumber (Double.NaN, Double.NaN);
		}
	}
	
	@Override
	double value(String [] var, double [] val) {
		// The value is obtained by evaluating the left and right
		// operands and combining the values with the operator.
		double x = left.value(var, val);
		double y = right.value(var, val);
		
		switch(op) {
		case '+': return x + y;
		case '-': return x - y;
		case '*': return x * y;
		case '/': return x / y;
		case '^': return Math.pow(x, y);
		default: return Double.NaN;	// bad operator!
		}
	}

	@Override
	StringBuffer printPostfix() {
		// evaluate the operands
		StringBuffer x = left.printPostfix();
		StringBuffer y = right.printPostfix();
		return new StringBuffer (x + " " + y + " " + op);
	}

	@Override
	StringBuffer printInfix() {
		String x = left.printInfix().toString();
		String y = right.printInfix().toString();
		
		switch (op) {
		case '+': case '-':
			String str = x + " " + op + " " + y;
			if (str.startsWith("0.0 + "))
				str = str.substring(6, str.length());
			if (str.startsWith(" + "))
				str = str.substring(3, str.length());
			if (str.endsWith(" + ") || str.endsWith(" - "))
				str = str.substring(0, str.length() - 3);
			if (str.endsWith(" + 0.0") || str.endsWith(" - 0.0"))
				str = str.substring(0, str.length() - 6);
			return new StringBuffer (str);
		case '*': case '/':
			if (x.equals("0.0") || y.equals("0.0") || 
					x.equals("") || y.equals(""))
				return new StringBuffer ("0.0");
			if (x.equals("1.0") && (op == '*'))
				return new StringBuffer (y);
			if (y.equals("1.0"))
				return new StringBuffer (x);
			return new StringBuffer (
				(hasOperator(x) ? " (" + x + ")" : x) +
				(hasOperator(y) ? 
						" " + op + " (" + y + ")" : " " + op + " " + y));
			
		case '%':
			if (x.length() != 0 && y.length() != 0)
				return new StringBuffer (
						(hasOperator(x) ? "(" + x + ") ": x + " ") + op +
						(hasOperator(y) ? " (" + y + ") ": " " + y));
			if (y.toString().isEmpty() && op == '^')
				return new StringBuffer ("1");
			if (y.toString().isEmpty() && op == '%')
				return new StringBuffer("0");
			throw new IllegalArgumentException ("Invalid operator");
		case '^':
			if (y.equals("0.0") || y.equals("") || x.equals("1.0"))
				return new StringBuffer ("1.0");
			if (y.equals("1.0"))
				return new StringBuffer (x);
			if (x.equals("0.0"))
				return new StringBuffer ("0.0");
			return new StringBuffer (
					(hasOperator(x) ? "(" + x + ")": x) + op +
					(hasOperator(y) ? " (" + y + ")": " " + y));
		default:
			throw new IllegalArgumentException ("Invalid operator");
		}
	}

	@Override
	ExpNode derivative(String x) {
		// Apply the derivative formulas
		switch (op) {
		case '+': case '-':
			return BinOpNode.make (op, left.derivative(x), right.derivative(x));
		case '*':
			return BinOpNode.make ('+', 
					BinOpNode.make ('*', left, right.derivative(x)),
					BinOpNode.make ('*', right, left.derivative(x)));
		case '/':
			return BinOpNode.make ('/',
					BinOpNode.make ('-',
							BinOpNode.make ('*', right, left.derivative(x)),
							BinOpNode.make ('*', left, right.derivative(x))),
							BinOpNode.make ('*', right, right));
		case '^':
			// for x ^ n, return n * x ^ n - 1
			if (!right.wrt().equals(x))
				return BinOpNode.make ('*', right, BinOpNode.make ('^', left,
						BinOpNode.make ('-', right, new ConstantNode(1))));
			return BinOpNode.make ('*', right.derivative(x),
					BinOpNode.make ('*', this, new LnNode(left)));
		default:
			return null;
		}
	}
}

class LnNode extends ExpNode {
	// instance variable
	ExpNode operand;
	
	// constructor for type ln(x)
	LnNode (ExpNode x) {
		this.operand = x;
	}
	
	@Override boolean isZero () {return operand.isUnity(); }
	@Override boolean isUnity () { return false;}
	
	@Override
	double value(String [] var, double [] val) {
		// TODO Auto-generated method stub
		return Math.log(operand.value(var, val));
	}
	
	@Override
	ComplexNumber complexValue (String [] var, double [] val) {
		return new ComplexNumber (value(var, val));
	}
	
	@Override
	StringBuffer printPostfix() {
		// TODO Auto-generated method stub
		return new StringBuffer (operand.printPostfix() + " ln ");
	}

	@Override
	StringBuffer printInfix() {
		// TODO Auto-generated method stub
		return new StringBuffer ("ln (" + operand.printInfix() + ")");
	}

	@Override
	ExpNode derivative(String wrt) {
		// TODO Auto-generated method stub
		return BinOpNode.make ('/', operand.derivative(wrt), operand);
	}
}

class ExponentialNode extends ExpNode {
	// instance variable
	ExpNode operand;
	
	// constructor for type exp(fx)
	ExponentialNode (ExpNode fx) {
		this.operand = fx;
	}
	
	@Override boolean isZero () { return false; }
	@Override boolean isUnity () { return operand.isZero(); }
	
	@Override
	double value(String [] var, double [] val) {
		// TODO Auto-generated method stub
		return Math.exp(operand.value(var, val));
	}
	
	@Override
	ComplexNumber complexValue (String [] var, double [] val) {
		ComplexNumber base = operand.complexValue(var, val);
		double im = base.getImaginaryNumber();
		ComplexNumber c = new ComplexNumber (Math.cos(im), Math.sin(im));
		return c.product(base.getRealNumber());
	}

	@Override
	StringBuffer printPostfix() {
		// TODO Auto-generated method stub
		return new StringBuffer (operand.printPostfix() + " exp ");
	}

	@Override
	StringBuffer printInfix() {
		// TODO Auto-generated method stub
		if (operand.printInfix().toString().equals("0.0"))
			return new StringBuffer ("1.0");
		return new StringBuffer ("exp (" + operand.printInfix() + ")");
	}

	@Override
	ExpNode derivative(String x) {
		// TODO Auto-generated method stub
		return BinOpNode.make ('*', operand.derivative(x), this);
	}
}

class SinNode extends ExpNode {
	// instance variables
	ExpNode operand;
	
	// constructor of type sin(fx)
	SinNode (ExpNode fx) {
		this.operand = fx;
	}
	
	@Override boolean isZero () { return operand.isZero(); }
	@Override boolean isUnity () { return false; }
	
	@Override
	double value(String [] var, double [] val) {
		// TODO Auto-generated method stub
		return Math.sin(operand.value(var, val));
	}
	
	@Override
	ComplexNumber complexValue (String [] var, double [] val) {
		return new ComplexNumber (value(var, val));
	}

	@Override
	StringBuffer printPostfix() {
		// TODO Auto-generated method stub
		return new StringBuffer (operand.printPostfix() + " sin");
	}

	@Override
	StringBuffer printInfix() {
		// TODO Auto-generated method stub
		return new StringBuffer ("sin (" + operand.printInfix() + ")");
	}

	@Override
	ExpNode derivative(String x) {
		// TODO Auto-generated method stub
		return BinOpNode.make ('*', operand.derivative(x), new CosNode(operand));
	}
}

class CosNode extends ExpNode {
	// instance variables
	ExpNode operand;
	
	// constructor of type sin(fx)
	CosNode (ExpNode fx) {
		this.operand = fx;
	}
	
	@Override boolean isZero () { return false; }
	@Override boolean isUnity () { return operand.isZero(); }
	
	@Override
	double value(String [] var, double [] val) {
		// TODO Auto-generated method stub
		return Math.cos(operand.value(var, val));
	}
	
	@Override
	ComplexNumber complexValue (String [] var, double [] val) {
		return new ComplexNumber (value(var, val));
	}

	@Override
	StringBuffer printPostfix() {
		// TODO Auto-generated method stub
		return new StringBuffer (operand.printPostfix() + " cos");
	}

	@Override
	StringBuffer printInfix() {
		// TODO Auto-generated method stub
		return new StringBuffer ("cos (" + operand.printInfix() + ")");
	}

	@Override
	ExpNode derivative(String x) {
		// TODO Auto-generated method stub
		return BinOpNode.make ('*', operand.derivative(x), 
				OperandNode.make (-1.0, new SinNode(operand), 1.0));
	}
}

class TanNode extends ExpNode {
	// instance variables
	ExpNode operand;
	
	// constructor of type sin(fx)
	TanNode (ExpNode fx) {
		this.operand = fx;
	}
	
	@Override boolean isZero () { return operand.isZero(); }
	@Override boolean isUnity () { return false; }
	
	@Override
	double value(String [] var, double [] val) {
		// TODO Auto-generated method stub
		return Math.tan(operand.value(var, val));
	}
	
	@Override
	ComplexNumber complexValue (String [] var, double [] val) {
		return new ComplexNumber (value(var, val));
	}

	@Override
	StringBuffer printPostfix() {
		// TODO Auto-generated method stub
		return new StringBuffer (operand.printPostfix() + " tan");
	}

	@Override
	StringBuffer printInfix() {
		// TODO Auto-generated method stub
		return new StringBuffer ("tan (" + operand.printInfix() + ")");
	}

	@Override
	ExpNode derivative(String x) {
		// TODO Auto-generated method stub
		return BinOpNode.make ('*', operand.derivative(x),
				OperandNode.make (1.0, new SecNode(operand), 2));
	}
}

class CosecNode extends ExpNode {
	// instance variables
	ExpNode operand;
	
	// constructor of type sin(fx)
	CosecNode (ExpNode fx) {
		this.operand = fx;
	}
	
	@Override boolean isZero () { return false; }
	@Override boolean isUnity () { return operand.isUnity(); }
	
	@Override
	double value(String [] var, double [] val) {
		// TODO Auto-generated method stub
		return Math.pow(Math.sin(operand.value(var, val)), -1.0);
	}
	
	@Override
	ComplexNumber complexValue (String [] var, double [] val) {
		return new ComplexNumber (value(var, val));
	}

	@Override
	StringBuffer printPostfix() {
		// TODO Auto-generated method stub
		return new StringBuffer (operand.printPostfix() + " cosec");
	}

	@Override
	StringBuffer printInfix() {
		// TODO Auto-generated method stub
		return new StringBuffer ("cosec (" + operand.printInfix() + ")");
	}

	@Override
	ExpNode derivative(String x) {
		// TODO Auto-generated method stub
		return BinOpNode.make ('*', operand.derivative(x), 
				OperandNode.make (-1.0,	BinOpNode.make ('*', 
						new CotNode (operand), this), 1.0));
	}
}

class SecNode extends ExpNode {
	// instance variables
	ExpNode operand;
	
	// constructor of type sin(fx)
	SecNode (ExpNode fx) {
		this.operand = fx;
	}
	
	@Override boolean isZero () { return false; }
	@Override boolean isUnity () { return operand.isUnity(); }
	
	@Override
	double value(String [] var, double [] val) {
		// TODO Auto-generated method stub
		return Math.pow(Math.cos(operand.value(var, val)), -1.0);
	}
	
	@Override
	ComplexNumber complexValue (String [] var, double [] val) {
		return new ComplexNumber (value(var, val));
	}

	@Override
	StringBuffer printPostfix() {
		// TODO Auto-generated method stub
		return new StringBuffer (operand.printPostfix() + " sec");
	}

	@Override
	StringBuffer printInfix() {
		// TODO Auto-generated method stub
		return new StringBuffer ("sec (" + operand.printInfix() + ")");
	}

	@Override
	ExpNode derivative(String x) {
		// TODO Auto-generated method stub
		return BinOpNode.make ('*', operand.derivative(x),
				BinOpNode.make ('*', new TanNode (operand), this));
	}
}

class CotNode extends ExpNode {
	// instance variables
	ExpNode operand;
	
	// constructor of type sin(fx)
	CotNode (ExpNode fx) {
		this.operand = fx;
	}
	
	@Override boolean isZero () { return false; }
	@Override boolean isUnity () { return operand.isUnity(); }
	
	@Override
	double value(String [] var, double [] val) {
		// TODO Auto-generated method stub
		return Math.pow(Math.tan(operand.value(var, val)), -1.0);
	}
	
	@Override
	ComplexNumber complexValue (String [] var, double [] val) {
		return new ComplexNumber (value(var, val));
	}

	@Override
	StringBuffer printPostfix() {
		// TODO Auto-generated method stub
		return new StringBuffer (operand.printPostfix() + " cot");
	}

	@Override
	StringBuffer printInfix() {
		// TODO Auto-generated method stub
		return new StringBuffer ("cot (" + operand.printInfix() + ")");
	}

	@Override
	ExpNode derivative(String x) {
		// TODO Auto-generated method stub
		return BinOpNode.make ('*', operand.derivative(x), 
				OperandNode.make (-1.0, new CosecNode (operand), 2.0));
	}
}

class SinhNode extends ExpNode {
	// instance variables
	ExpNode operand;
	
	// constructor of type sin(fx)
	SinhNode (ExpNode fx) {
		this.operand = fx;
	}
	
	@Override boolean isZero () { return operand.isZero(); }
	@Override boolean isUnity () { return false; }
	
	@Override
	double value(String [] var, double [] val) {
		// TODO Auto-generated method stub
		return Math.sinh(operand.value(var, val));
	}
	
	@Override
	ComplexNumber complexValue (String [] var, double [] val) {
		return new ComplexNumber (value(var, val));
	}

	@Override
	StringBuffer printPostfix() {
		// TODO Auto-generated method stub
		return new StringBuffer (operand.printPostfix() + " sinh");
	}

	@Override
	StringBuffer printInfix() {
		// TODO Auto-generated method stub
		return new StringBuffer ("sinh (" + operand.printInfix() + ")");
	}

	@Override
	ExpNode derivative(String x) {
		// TODO Auto-generated method stub
		return BinOpNode.make ('*', operand.derivative(x), 
				new CoshNode (operand));
	}
}
class CoshNode extends ExpNode {
	// instance variables
	ExpNode operand;
	
	// constructor of type sin(fx)
	CoshNode (ExpNode fx) {
		this.operand = fx;
	}
	
	@Override boolean isZero () { return false; }
	@Override boolean isUnity () { return operand.isZero(); }
	
	@Override
	double value(String [] var, double [] val) {
		// TODO Auto-generated method stub
		return Math.cosh(operand.value(var, val));
	}
	
	@Override
	ComplexNumber complexValue (String [] var, double [] val) {
		return new ComplexNumber (value(var, val));
	}

	@Override
	StringBuffer printPostfix() {
		// TODO Auto-generated method stub
		return new StringBuffer (operand.printPostfix() + " cosh");
	}

	@Override
	StringBuffer printInfix() {
		// TODO Auto-generated method stub
		return new StringBuffer ("cosh (" + operand.printInfix() + ")");
	}

	@Override
	ExpNode derivative(String x) {
		// TODO Auto-generated method stub
		return BinOpNode.make ('*', operand.derivative(x),
				new SinhNode (operand));
	}
}
class TanhNode extends ExpNode {
	// instance variables
	ExpNode operand;
	
	// constructor of type sin(fx)
	TanhNode (ExpNode fx) {
		this.operand = fx;
	}
	
	@Override boolean isZero () { return operand.isZero(); }
	@Override boolean isUnity () { return false; }
	
	@Override
	double value(String [] var, double [] val) {
		// TODO Auto-generated method stub
		return Math.tanh(operand.value(var, val));
	}
	
	@Override
	ComplexNumber complexValue (String [] var, double [] val) {
		return new ComplexNumber (value(var, val));
	}

	@Override
	StringBuffer printPostfix() {
		// TODO Auto-generated method stub
		return new StringBuffer (operand.printPostfix() + " tanh");
	}

	@Override
	StringBuffer printInfix() {
		// TODO Auto-generated method stub
		return new StringBuffer ("tanh (" + operand.printInfix() + ")");
	}

	@Override
	ExpNode derivative(String x) {
		// TODO Auto-generated method stub
		return BinOpNode.make ('*', operand.derivative(x),
				OperandNode.make (1.0, new SechNode (operand), 2.0));
	}
}

class SechNode extends ExpNode {
	// instance variables
	ExpNode operand;
	
	// constructor of type sin(fx)
	SechNode (ExpNode fx) {
		this.operand = fx;
	}
	
	@Override boolean isZero () { return false; }
	@Override boolean isUnity () { return operand.isUnity(); }
	
	@Override
	double value(String [] var, double [] val) {
		// TODO Auto-generated method stub
		return Math.pow(Math.cosh(operand.value(var, val)), -1.0);
	}
	
	@Override
	ComplexNumber complexValue (String [] var, double [] val) {
		return new ComplexNumber (value(var, val));
	}

	@Override
	StringBuffer printPostfix() {
		// TODO Auto-generated method stub
		return new StringBuffer (operand.printPostfix() + " sech");
	}

	@Override
	StringBuffer printInfix() {
		// TODO Auto-generated method stub
		return new StringBuffer ("sech (" + operand.printInfix() + ")");
	}

	@Override
	ExpNode derivative(String x) {
		// TODO Auto-generated method stub
		return BinOpNode.make ('*', operand.derivative(x), 
				OperandNode.make (-1.0, BinOpNode.make ('*', 
						new TanhNode (operand), this), 1.0));
	}
}

class CosechNode extends ExpNode {
	// instance variables
	ExpNode operand;
	
	// constructor of type sin(fx)
	CosechNode (ExpNode fx) {
		this.operand = fx;
	}
	
	@Override boolean isZero () { return false; }
	@Override boolean isUnity () { return operand.isUnity(); }
	
	@Override
	double value(String [] var, double [] val) {
		// TODO Auto-generated method stub
		return Math.pow(Math.sinh(operand.value(var, val)), -1.0);
	}
	
	@Override
	ComplexNumber complexValue (String [] var, double [] val) {
		return new ComplexNumber (value(var, val));
	}

	@Override
	StringBuffer printPostfix() {
		// TODO Auto-generated method stub
		return new StringBuffer (operand.printPostfix() + " cosec");
	}

	@Override
	StringBuffer printInfix() {
		// TODO Auto-generated method stub
		return new StringBuffer ("cosec (" + operand.printInfix() + ")");
	}

	@Override
	ExpNode derivative(String x) {
		// TODO Auto-generated method stub
		return BinOpNode.make ('*', operand.derivative(x), 
				OperandNode.make (-1.0, BinOpNode.make ('*', 
						new CothNode (operand), this), 1.0));
	}
}
class CothNode extends ExpNode {
	// instance variables
	ExpNode operand;
	
	// constructor of type sin(fx)
	CothNode (ExpNode fx) {
		this.operand = fx;
	}
	
	@Override boolean isZero () { return false; }
	@Override boolean isUnity () { return operand.isUnity(); }
	
	@Override
	double value(String [] var, double [] val) {
		// TODO Auto-generated method stub
		return Math.pow(Math.tanh(operand.value(var, val)), -1.0);
	}
	
	@Override
	ComplexNumber complexValue (String [] var, double [] val) {
		return new ComplexNumber (value(var, val));
	}

	@Override
	StringBuffer printPostfix() {
		// TODO Auto-generated method stub
		return new StringBuffer (operand.printPostfix() + " coth");
	}

	@Override
	StringBuffer printInfix() {
		// TODO Auto-generated method stub
		return new StringBuffer ("coth (" + operand.printInfix() + ")");
	}

	@Override
	ExpNode derivative(String x) {
		// TODO Auto-generated method stub
		return BinOpNode.make ('*', operand.derivative(x),
				OperandNode.make (-1.0, new CosechNode(operand), 2.0));
	}
}