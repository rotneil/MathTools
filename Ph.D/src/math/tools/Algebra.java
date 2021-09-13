package math.tools;

/**
 * This class contains static methods which are used as utility
 * methods for performing operations on polynomial expressions and
 * equations.
 * 
 * @author Oluwafemi Nehemiah
 * rotneil@yahoo.com
 *
 */
public final class Algebra {
	
	/**
	 * This method uses the Runge Kutta Fourth order numerical integration
	 * scheme to return the next Poincare Section crossing of a function 
	 * expressed as ExpNode fx. The accuracy of the crossing is determined 
	 * by err which is evaluated by using the halving method (i.e 
	 * Bisection method) for the poicare section interpolation.
	 * @param fx The function equation
	 * @param xIndex The indices of the independent variables in String array
	 * variables
	 * @param sectionIndex The index of the Poincare Section variable q(x) 
	 * in var
	 * @param var The variables that are defined in the state function fx
	 * @param val The corresponding variable values
	 * @param The index of the time signature in val
	 * @param dt The step interval used in the Runge Kutta integration scheme
	 * @param err This is the maximum allowable error in the evaluation of the
	 * Poincare section crossing
	 * 
	 * @return The Poicare Crossing coordinates with the time for that 
	 * crossing (to the accuracy of err).
	 * 
	 * @throws IllegalStateException iff the orbit goes into an equilibrium
	 * point during the Poincare crossing scheme
	 */
	public static double [] getNextCrossing (ExpNode [] fx, int [] xIndex,
			int sectionIndex, String [] var, double [] val, 
			int timeIndex, double dt, double err)
	{
		// local variables
		double [] v = rungeKutta (fx, xIndex, var, val, dt);
		double [] vv = rungeKutta (fx, xIndex, var, v, dt);
		int count = 1;
		while (!(v[xIndex[sectionIndex]] <= 0.0 && vv[xIndex[sectionIndex]] > 0.0)){
			// check that the orbit is not at equilibrium point
			if (Math.abs(v[xIndex[0]] - vv[xIndex[0]]) < 1.0E-12 ||
					Double.isNaN(vv[sectionIndex]))
				throw new IllegalStateException (
						"The orbit is at equilibrium point.");
			v = vv;
			vv = rungeKutta (fx, xIndex, var, v, dt);
			
			++count;
		}
		
		// set the rough time estimate
		if (timeIndex != -1)
			v[timeIndex] = dt * count;
		
		// check if orbit is exactly on the section
		if (v[xIndex[sectionIndex]] == 0.0)
			return v;
		
		return timeHalf (fx, xIndex, sectionIndex, var, v, timeIndex, dt, err);
	}
	
	// timeHalf method for Poincare section mapping
	public static double [] timeHalf (ExpNode [] fx, int [] xIndex, int sectionIndex, 
			String [] var, double [] val, int timeIndex, double dt, double err)
	{
		// local variable
		int count = 0;
		double tt = 0.0;
		double dtt = dt * 0.5;
		double [] vv = rungeKutta (fx, xIndex, var, val, dtt);
		while (Math.abs(vv[xIndex[sectionIndex]]) > err && ++count < 200) {
			if (vv[xIndex[sectionIndex]] > 0.0)
				dt = dtt;
			else
				tt = dtt;
			
			dtt = (dt + tt) * 0.5;
			vv = rungeKutta (fx, xIndex, var, val, dtt);
		}
		
		// get the time index and put t into value
		if (timeIndex != -1)
			vv[timeIndex] += dtt;
		
		return vv;
	}
	
	/**
	 * This is a standard RungeKutta method for integration of 
	 * an array of equation as ExpNode. An instance variation is computed
	 * and returned
	 * 
	 * @param fx The equations to integrate
	 * @param xIndex The indices of the independent variables in the array 
	 * of variables var
	 * @param var The variable that are defined in the ExpNode fx
	 * @param val The corresponding values of the variables defined in var
	 * @param dt The integration interval for the rungekutta scheme.
	 * @return The double []
	 */
	public static double [] rungeKutta (ExpNode [] fx, int [] xIndex, 
			String [] var, double [] val, double dt)
 	{
 		int l = fx.length;
 		double [] c1, c2, c3, c4, temp;
 		
 		// initialize the intermediate steps
 		c1 = new double [l];
 		c2 = new double [l];
 		c3 = new double [l];
 		c4 = new double [l];
 		temp = new double [val.length];
 		
 		// copy the content of val into the temporary array
 		for (int i = 0; i < temp.length; ++i)
 			temp[i] = val[i];
 		
 		// perform the intemediate integration steps
 		c1 = Matrix.values(fx, var, val);
 		
 		for (int i = 0; i < l; ++i)
 			temp[xIndex[i]] = val[xIndex[i]] + dt * c1[i] / 2;
 		c2 = Matrix.values(fx, var, temp);
 		
 		for (int i = 0; i < l; ++i)
 			temp[xIndex[i]] = val[xIndex[i]] + dt * c2[i] / 2;
 		c3 = Matrix.values(fx, var, temp);
 		
 		for (int i = 0; i < l; ++i)
 			temp[xIndex[i]] = val[xIndex[i]] + dt * c3[i];
 		c4 = Matrix.values(fx, var, temp);
 		
 		for (int i = 0; i < l; ++i)
 			temp[xIndex[i]] = val[xIndex[i]] + 
 				dt * (c1[i] + 2 * (c2[i] + c3[i]) + c4[i]) / 6.0;
 		
 		return temp;
 	}
 	
	/**
	 * This method uses the native RungeKutta Fourth order integration scheme
	 * to calculate the variables of the next Poincare Section crossing
	 * 
	 * @param section The coordinate about which the Poincare Section is determined
	 * @param v The initial values used in the integration process whose index 
	 * consists of the time index as the last.
	 * 
	 * @param p The parameters of the state function
	 * @param dt The intermediate step of the RungeKutta scheme
	 * @param err The error estimate for the Poincare crossing
	 * 
	 * @return
	 */
 	public static double [] getNextCrossing (int section, 
 			double [] v, double [] p, double dt, double err)
	{
		// local variables
		double [] v1 = rungeKutta (v, p, dt);
		double [] v2 = rungeKutta (v1, p, dt);
		int count = 1;
		while (!(v1[1] < 0.0 && v2[1] >= 0.0)){
			// check that the orbit is not at equilibrium point
			if (Math.abs(v1[1] - v2[1]) < 1.0E-12)
				throw new IllegalStateException (
						"The orbit is at equilibrium point.");
			v1 = v2;
			v2 = rungeKutta (v1, p, dt);
			
			++count;
		}
		
		// set the rough time estimate
		v1[v1.length - 1] = count * dt;
		
		return timeHalf (section, v1, p, dt, err);
	}
	
 	/**
 	 * This is method is based on bisection method for refinning the values of
 	 * the state variables at Poincare Section Crossing.
 	 * @param section The coordinate about which the Poincare Section is taken
 	 * @param v
 	 * @param p
 	 * @param dt
 	 * @param err
 	 * 
 	 * @return
 	 */
	// timeHalf method for Poincare section mapping
	public static double [] timeHalf (
			int section, double [] v, double [] p, double dt, double err)
	{
		// local variable
		int count = 0;
		double tt = 0.0;
		double dtt = dt * 0.5;
		double [] vv = rungeKutta (v, p, dtt);
		while (Math.abs(vv[section]) > err && ++count < 200) {
			if (vv[section] > 0.0)
				dt = dtt;
			else
				tt = dtt;
			dtt = (dt + tt) * 0.5;
			vv = rungeKutta (v, p, dtt);
		}
		
		// get the time index and put t into value
		vv[vv.length - 1] += dtt;
		
		return vv;
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
	public static double [] rungeKutta (double [] v, double [] p, double dt)
 	{
 		int l = v.length;	// putting in mind the time index
 		double [] c1, c2, c3, c4;
 		
 		// initialize the intermediate steps
 		c1 = new double [l];
 		c2 = new double [l];
 		c3 = new double [l];
 		c4 = new double [l];
 		
 		c1 = f(v, p[0], p[1]);
 		
 		for (int i = 0; i < l; ++i)
 			c2[i] = v[i] + dt * c1[i] / 2;
 		c2 = f(c2, p[0], p[1]);
 		
 		for (int i = 0; i < l; ++i)
 			c3[i] = v[i] + dt * c2[i] / 2;
 		c3 = f(c3, p[0], p[1]);
 		
 		for (int i = 0; i < l; ++i)
 			c4[i] = v[i] + dt * c3[i];
 		c4 = f(c4, p[0], p[1]);
 		
 		for (int i = 0; i < l; ++i)
 			c1[i] = v[i] + dt * (c1[i] + 2 * (c2[i] + c3[i]) + c4[i]) / 6.0;
 		
 		return c1;
 	}
	
	/**
	 * This is the state variable definition which may include those of the
	 * first and/or second variation equation.
	 * @param v
	 * @param g
	 * @param k
	 * @return
	 */
	public static double [] f (double [] v, double g, double k)
	{
		// declare local variable
		double [] vv = new double [v.length];
		double sech = 1.0 / Math.cosh (g * v[0]);
		double tanh = Math.tanh (g * v[0]);
		
		vv[0] = -v[1] + tanh;
		vv[1] = v[0] - k * v[1];
		
		// variation
		vv[2] = g * sq(sech) * v[2] - v[4];
		vv[3] = g * sq(sech) * v[3] - v[5];
		vv[4] = v[2] - k * v[4];
		vv[5] = v[3] - k * v[5];
		
		vv[6] = g * sq(sech) * v[6] - v[7];
		vv[7] = v[6] - k * v[7] - v[1];
		
		vv[8] = g * sq(sech) * v[8] - 2 * tanh * sq(g * sech * v[2]) - v[11];
		vv[9] = g * sq(sech) * v[9] - 2 * tanh * sq(g * sech) * v[2] * v[3] - v[12];
		vv[10] = g * sq(sech) * v[10] - 2 * sq(g * sech * v[3]) * tanh - v[13];
		vv[11] = v[8] - k * v[11];
		vv[12] = v[9] - k * v[12];
		vv[13] = v[10] - k * v[13];
		
		vv[14] = g * sq(sech) * v[14] - 2 * sq(g * sech) * tanh * v[2] * v[6] - v[16];
		vv[15] = g * sq(sech) * v[15] - 2 * sq(g * sech) * tanh * v[3] * v[6] - v[17];
		vv[16] = v[14] - k * v[16] - v[4];
		vv[17] = v[15] - k * v[17] - v[5];
				
		return vv;
	}
 	
 	public static double sq (double x) {return x * x; }
	
	/**
	 * This method reads an infix expression and returns it 
	 * as a postfix expression. The expression may be composed of digits alone
	 * or an algebra expression.
	 * @param infix The arithmetics expression
	 * @return
	 */
	public static StringBuffer infixToPostfix (StringBuffer exp)
	{
		StringBuffer postfix = new StringBuffer();
		Stack <String> stack = new Stack <String> ();
		
		// push a left parenthesis into the stack and append ')' into infix
		stack.push("(");
		StringBuffer infix = new StringBuffer (exp);
		infix.append(')');
		
		// perform the recursive operation
		int index = 0;	// this is the read index in infix
		char current;
		while (!stack.isEmpty()) {
			// read the current character in infix
			current = infix.charAt(index++);
			
			// eliminate white spaces
			while (Character.isWhitespace(current))
				current = infix.charAt(index++);
			
			// if current is a digit
			if (Character.isDigit(current)) 
			{
				// read a number more than a unit
				StringBuffer number = new StringBuffer();
				number.append(current);
				while (Character.isDigit(infix.charAt(index)) || 
						infix.charAt(index) == '.')
					number.append(infix.charAt(index++));
				
				// check for concatenated digits-alphabet like 3x
				if (Character.isAlphabetic(infix.charAt(index))) {
					StringBuffer variable = new StringBuffer();
					variable.append(infix.charAt(index++));
					while (Character.isAlphabetic(infix.charAt(index)) || 
							Character.isDigit(infix.charAt(index)) || 
							infix.charAt(index) == '_')
						variable.append(infix.charAt(index++));
					
					// check if the variable is written with power char
					while (infix.charAt(index) == ' ') ++index;
					if (infix.charAt(index) == '^') {
						// increment infix index
						index++;
						
						// define the index
						StringBuffer numIndex = new StringBuffer ();
						
						// remove the white spaces
						while (infix.charAt(index) == ' ') index++;
						
						// read the next digits or variables
						if (Character.isDigit(infix.charAt(index)) ||
								Character.isAlphabetic(infix.charAt(index)) ||
								infix.charAt(index) == '.') 
						{
							numIndex.append(infix.charAt(index++));
							while (Character.isDigit(infix.charAt(index)) ||
								Character.isAlphabetic(infix.charAt(index)) ||
								infix.charAt(index) == '.')
								numIndex.append(infix.charAt(index++));
							
							// form the postfix expression
							postfix.append(number + " " + variable + " " +
									numIndex + " ^ * ");
						}
					} else {
						// if the digit is writen with a variable,
						// the digit and variable with a times to postfix
						postfix.append(number + " " + variable + " * ");
					}
				} else
					// append to postfix
					postfix.append(number + " ");
			}
			
			// if the current is for algebraic variable
			else if (Character.isAlphabetic(current)) 
			{
				// read a number more than a unit
				StringBuffer variable = new StringBuffer();
				variable.append(current);
				while (Character.isAlphabetic(infix.charAt(index)) || 
						Character.isDigit(infix.charAt(index)) ||
						infix.charAt(index) == '_')
					variable.append(infix.charAt(index++));
				
				// check whether the alphabet is a function
				if (isFunction (variable.toString()))
					// push current operator into stack
					stack.push(variable.toString());
				else {
					// check if the variable is written with power char
					while (infix.charAt(index) == ' ') ++index;
					if (infix.charAt(index) == '^') {
						// increment infix index
						index++;
						
						// define the index
						StringBuffer numIndex = new StringBuffer ();
						
						// remove the white spaces
						while (infix.charAt(index) == ' ') index++;
						
						// read the next digits or variables
						if (Character.isDigit(infix.charAt(index)) ||
								Character.isAlphabetic(infix.charAt(index)) ||
								infix.charAt(index) == '.') {
							numIndex.append(infix.charAt(index++));
							
							while (Character.isDigit(infix.charAt(index)) ||
								Character.isAlphabetic(infix.charAt(index)) ||
								infix.charAt(index) == '.')
								numIndex.append(infix.charAt(index++));
							
							// form the postfix expression
							postfix.append(variable + " " +	numIndex + " ^ ");
						}
					} else {
						// append variables to postfix
						postfix.append(variable + " ");
					}
				}
			}
			
			// current is right parathensis
			else if (current == '(')
				// push it into stack
				stack.push("(");
			else if (isOperator (current)) {
				// check for negative numbers
                if (current == '-' && Character.isDigit(infix.charAt(index + 1))) {
                    {
                        // read a number more than a unit
                        StringBuffer number = new StringBuffer();
                        number.append(current);
                        
                        // continue with the rest of the digits
                        while (Character.isDigit(infix.charAt(index)) || 
                                    infix.charAt(index) == '.')
                            number.append(infix.charAt(index++));

                        // check for concatenated digits-alphabet like 3x
                        if (Character.isAlphabetic(infix.charAt(index))) {
                            StringBuffer variable = new StringBuffer();
                            variable.append(infix.charAt(index++));
                            while (Character.isAlphabetic(infix.charAt(index)) || 
                                    Character.isDigit(infix.charAt(index)) || 
                                    infix.charAt(index) == '_')
                                variable.append(infix.charAt(index++));

                            // check if the variable is written with power char
                            while (infix.charAt(index) == ' ') ++index;
                            if (infix.charAt(index) == '^') {
                                // increment infix index
                                index++;

                                // define the index
                                StringBuffer numIndex = new StringBuffer ();

                                // remove the white spaces
                                while (infix.charAt(index) == ' ') index++;

                                // read the next digits or variables
                                if (Character.isDigit(infix.charAt(index)) ||
                                        Character.isAlphabetic(infix.charAt(index)) ||
                                        infix.charAt(index) == '.') 
                                {
                                    numIndex.append(infix.charAt(index++));
                                    while (Character.isDigit(infix.charAt(index)) ||
                                            Character.isAlphabetic(infix.charAt(index)) ||
                                            infix.charAt(index) == '.')
                                        numIndex.append(infix.charAt(index++));

                                    // form the postfix expression
                                    String arg = number + " " + variable + " " +
                                                    numIndex + " ^ * ";
                                    postfix.append(arg);
                                }
                            } else {
                                // if the digit is writen with a variable,
                                // the digit and variable with a times to postfix
                                postfix.append(number).append(" ").append(variable).append(" * ");
                            }
                        } else
                            // append to postfix
                            postfix.append(number).append(" ");
                    }
                } 
                
                // now enter normal operator input
                else {
                    // pop operators of higher or equal precedence from stack
                    while ((isOperator (stack.peek()) || isFunction(stack.peek()))
                                && precedence(current, stack.peek()))
                        postfix.append(stack.pop()).append(" ");
	                
	                // push current operator into stack
					stack.push("" + current);
                }
			}
			
			else if (current == ')') {
				// pop operators until a right parenthesis is at the top
				while ((isOperator (stack.peek()) || isFunction (stack.peek()))
						&& !stack.peek().equals("("))
					postfix.append(stack.pop() + " ");
				
				// pop and discard the last 
				stack.pop();
			}
		}
		return postfix;
	}
	

	
	/**
	 * This method converts an algebraic infix expression to an ExpressionNode.
	 * 
	 * @param infix The algebra expression which may or may not contain the
	 * independent variable defined by character x.
	 * @param variables The array that entails all the variables the expression 
	 * infix is consist of
	 * @param wrt That is derivative with-respect-to wrt
	 * @param values An array of values corresponding the the defined variables
	 * @return A Binary operation node
	 */
	public static ExpNode infixToBinOpNode (StringBuffer infix,
			String [] variables)
	{
		// declare the local variables
		StringBuffer postfix = infixToPostfix (infix);
		Stack <ExpNode> stack = new Stack <ExpNode> ();
		int index = 0;
		char current = postfix.charAt(index);
		ExpNode x, y;	// the operands
		
		// append right parathensis to postfix
		postfix.append(')');
		
		// look up the postfix expression
		while (current != ')') {
			// if current is a digit
			if (Character.isDigit(current)) {
				// read a number more than a unit
				StringBuffer number = new StringBuffer();
				number.append(current);
				while (Character.isDigit(postfix.charAt(index + 1)) || 
						postfix.charAt(index + 1) == '.')
					number.append(postfix.charAt(++index));
			
				// convert number to a constantNode and add to stack
				stack.push(
						new ConstantNode(Double.parseDouble(number.toString())));
			}
			
			// if the current is for algebraic variable
			else if (Character.isAlphabetic(current)) {
				// read a number more than a unit
				StringBuffer variable = new StringBuffer();
				variable.append(current);
				while (Character.isAlphabetic(postfix.charAt(index + 1)) || 
						Character.isDigit(postfix.charAt(index + 1)) || 
						postfix.charAt(index + 1) == '_')
					variable.append(postfix.charAt(++index));
				
				// push the variable to stack
				if (isDefined (variables, variable.toString()))	// option with the variable of x
					stack.push(new VariableNode (variable.toString()));
				// option for if the variable is a function
				else if (isFunction(variable.toString()))
					stack.push(assignFunctionNode (
							variable.toString(), stack.pop()));
				// option for other variables are considered as constant nodes
				else 
					throw new IllegalArgumentException (
							"Undefined variable " + variable.toString());
			}
			
			// if the current char is an operator
			else if (isOperator (current)) {
				// check for the case of a unary minus node
				if (stack.size() == 1) {
					if (current == '+')
						stack.push(new BinOpNode('+', new ConstantNode(0),
								stack.pop()));
					else if (current == '-')	// case '-'
						// stack the unary node
						stack.push(new OperandNode (-1.0, stack.pop(), 1.0));
				} else {
					// pop the next two values in the stack as the operands
					y = stack.pop();	// FILO
					x = stack.pop();
					
					// check for multiplication operator
					if (current == '*') {
						if (x instanceof ConstantNode && 
								y instanceof ConstantNode)
							stack.push(new ConstantNode (
									((ConstantNode) x).value() * 
									((ConstantNode) y).value()));
						else if (x instanceof ConstantNode && 
								y instanceof VariableNode)
							stack.push(new OperandNode (
									((ConstantNode) x).value(), y, 1.0));
						else if (x instanceof VariableNode && 
								y instanceof ConstantNode)
							stack.push(new OperandNode (
									((ConstantNode) y).value(), x, 1.0));
						else if (x instanceof VariableNode && 
								y instanceof VariableNode && 
								x.wrt().equals(y.wrt()) && !x.wrt().equals(""))
							stack.push(new OperandNode (1.0, x, 2.0));
						else if (x instanceof ConstantNode && 
								y instanceof OperandNode) {
							OperandNode yNode = (OperandNode) y;
							yNode.coef *= ((ConstantNode) x).value();
							stack.push(yNode);
						} else if (x instanceof OperandNode && 
								y instanceof ConstantNode) {
							OperandNode xNode = (OperandNode) x;
							xNode.coef *= ((ConstantNode) y).value();
							stack.push(xNode);
						} else if (x instanceof VariableNode && 
								y instanceof VariableNode && 
								x.wrt().equals(y.wrt()) && !x.wrt().equals(""))
							stack.push(new OperandNode (1.0, x, 2.0));
						else if (x instanceof VariableNode && 
								y instanceof OperandNode && 
								x.wrt().equals(y.wrt()) && !x.wrt().equals(""))
							stack.push(new OperandNode (
									((OperandNode) y).coef, x, 
									((OperandNode) y).index + 1.0));
						else if (x instanceof OperandNode && 
								y instanceof VariableNode && 
								x.wrt().equals(y.wrt()) && !x.wrt().equals(""))
							stack.push(new OperandNode (
									((OperandNode) x).coef, y, 
									((OperandNode) x).index + 1.0));
						else if (x instanceof OperandNode && 
								y instanceof OperandNode && 
								x.wrt().equals(y.wrt()) && !x.wrt().equals("")) {
							OperandNode xNode = (OperandNode) x;
							OperandNode yNode = (OperandNode) y;
							xNode.coef *= yNode.coef;
							xNode.index += yNode.index;
							stack.push(xNode);
						} else
							stack.push(new BinOpNode(current, x, y));
					} else if (current == '^' && 
							x instanceof OperandNode && 
							y instanceof ConstantNode) {
						// downcast x
						OperandNode opX = (OperandNode) x;
						double value = ((ConstantNode) y).value();
						
						// stack a new OperandNode with the new coefficients,
						// and index
						stack.push(new OperandNode (
								Math.pow(opX.coef, value), opX.x, value));
					} else if (current == '^' && 
							x instanceof VariableNode && 
							y instanceof ConstantNode) 
						stack.push(new OperandNode (1.0, x, 
								((ConstantNode) y).value()));
					else					
						// stack an instance of BinOpNode with the defined operator
						stack.push(new BinOpNode (current, x, y));
				}
			}
			
			// read the next character
			current = postfix.charAt(++index);
			while (current == ' ')
				current = postfix.charAt(++index);
		}
		
		// check that stack is left with just one element
		if (stack.size() > 1)
			throw new IllegalArgumentException ("Invalid equation syntax!");
		
		return stack.pop();
	}
	
	/**
	 * This method computes a postfix that is composed entirely of digits and
	 * returns the value as a double. This method throws an uncheck exception
	 * IllegalArgumentException whenever a charactor is found. See 
	 * computePostfix (postfix, variable[], values[])
	 * @param postfix
	 * @return
	 */
	public static double computePostfix (StringBuffer postfix)
	{
		// declare the local variables
		Stack <Double> stack = new Stack <Double> ();
		int index = 0;
		char current = postfix.charAt(index);
		double x, y;	// the operands
		
		// append right parathensis to postfix
		postfix.append(')');
		
		// look up the postfix expression
		while (current != ')') {
			// if current is a digit
			if (Character.isDigit(current)) {
				// read a number more than a unit
				StringBuffer number = new StringBuffer();
				number.append(current);
				while (Character.isDigit(postfix.charAt(index + 1)) || 
						postfix.charAt(index + 1) == '.')
					number.append(postfix.charAt(++index));
				
				// convert number to double and add to stack
				double value = Double.parseDouble(number.toString());
				stack.push(value);
			}
			
			// if the current is for algebraic variable
			else if (Character.isAlphabetic(current)) {
				// check for the possibility of a function
				StringBuffer variable = new StringBuffer();
				variable.append(current);
				while (Character.isAlphabetic(postfix.charAt(index)) || 
						Character.isDigit(postfix.charAt(index)) || 
						postfix.charAt(index) == '_')
					variable.append(postfix.charAt(index++));
				
				if (isFunction(variable.toString()))
					stack.push(calculate (variable.toString(), stack.pop()));
				else
					// throw an exception and prompt to use correct method
					throw new IllegalArgumentException ("Invalid character found." +
						" Use the right method call for the operation");
			}
			
			// if the current char is an operator
			else if (isOperator (current)) {
				// check for negative numbers
                if (current == '-' && Character.isDigit(postfix.charAt(index))){
                    // read a number more than a unit
                    StringBuilder number = new StringBuilder();
                    number.append(current);
                    
                    while (Character.isDigit(postfix.charAt(index + 1)) || 
                                    postfix.charAt(index + 1) == '.')
                            number.append(postfix.charAt(++index));

                    // convert number to double and add to stack
                    double value = Double.parseDouble(number.toString());
                    stack.push(value);
                } else {
                    // pop the next two values in the stack as the operands
                    y = stack.pop();
                    x = stack.pop();

                    // compute their operation and stack the result
                    stack.push(calculate (current, x, y));
                }
			}
			
			// read the next character
			current = postfix.charAt(++index);
		}
		
		return stack.pop();
	}
	
	
	/**
	 * This method computes the value of an infix algebraic expression
	 * where the variables are defined in variables[] and their 
	 * corresponding values in value.
	 * @param algebra The infix expression which is first converted to
	 * postfix and then evaluated on a stack machine
	 * @param variables The String array that has all the variable defined in it.
	 * @param values This is the corresponding values of the variables in the
	 * infix algebraic expression.
	 * @return
	 */
	public static double computeAlgebra (StringBuffer infix, String [] variables,
			double [] values)
	{
		if (variables.length != values.length)
			throw new IllegalArgumentException ("The length of array var must be " +
					"equal to the length of double value");
		
		// covert infix expression to postfix
		StringBuffer postfix = infixToPostfix (infix);
		
		// declare the local variables
		Stack <Double> stack = new Stack <Double> ();
		int index = 0;
		char current = postfix.charAt(index);
		double x, y;	// the operands
		
		// append right parathensis to postfix
		postfix.append(')');
		
		// look up the postfix expression
		while (current != ')') {
			// if current is a digit
			if (Character.isDigit(current)) {
				// read a number more than a unit
				StringBuffer number = new StringBuffer();
				number.append(current);
				while (Character.isDigit(postfix.charAt(index + 1)) || 
						postfix.charAt(index + 1) == '.')
					number.append(postfix.charAt(++index));
				
				// convert number to double and add to stack
				double value = Double.parseDouble(number.toString());
				stack.push(value);
			}
			
			// if the current is for algebraic variable
			else if (Character.isAlphabetic(current)) {
				// read a number more than a unit
				StringBuffer variable = new StringBuffer();
				variable.append(current);
				while (Character.isAlphabetic(postfix.charAt(index + 1)) || 
						Character.isDigit(postfix.charAt(index + 1)) || 
						postfix.charAt(index + 1) == '_')
					variable.append(postfix.charAt(++index));
				
				// push the varible to stack
				if (isDefined(variables, variable.toString()))
				stack.push (
					ExpNode.getValue(variable.toString(), variables, values));
			}
			
			// if the current char is an operator
			else if (isOperator (current)) {
				// pop the next two values in the stack as the operands
				y = stack.pop();
				x = stack.pop();
				
				// compute their operation and stack the result
				stack.push(calculate (current, x, y));
			}
			
			// read the next character
			current = postfix.charAt(++index);
		}
		
		return stack.pop();
	}
	
	private static boolean isDefined (String [] variables, String wrt)
	{
		for (String var : variables)
			if (var.equals(wrt))
				return true;
		return false;
	}
	
	// utiility method that does the actual calculation of the operations
	private static double calculate (char operator, double x, double y)
	{
		// perform computation based on the operator
		switch (operator) {
		case '-':
			return (x - y);
		case '+':
			return (x + y);
		case '*':
			return (x * y);
		case '/':
			return (x / y);
		case '^':
			return Math.pow(x, y);
		case '%':
			return (x % y);
		default:
			throw new UnknownOperatorException ("'" + operator + "' is invalid");
		}
	}
	
	// method to perform function calculation
	private static double calculate (String fx, double x)
	{
		if (fx.equals("ln")) return Math.log(x);
		else if (fx.equals("exp")) return Math.exp(x);
		else if (fx.equals("sin")) return Math.sin(x);
		else if (fx.equals("cos")) return Math.cos(x);
		else if (fx.equals("tan")) return Math.tan(x);
		else if (fx.equals("sinh")) return Math.sinh(x);
		else if (fx.equals("cosh")) return Math.cosh(x);
		else if (fx.equals("tanh")) return Math.tanh(x);
		else if (fx.equals("cosec")) return Math.pow(Math.sin(x), -1.0);
		else if (fx.equals("sec")) return Math.pow(Math.cos(x), -1.0);
		else if (fx.equals("cot")) return Math.pow(Math.tan(x), -1.0);
		else if (fx.equals("cosech")) return Math.pow(Math.sinh(x), -1.0);
		else if (fx.equals("sech")) return Math.pow(Math.cosh(x), -1.0);
		else if (fx.equals("coth")) return Math.pow(Math.tanh(x), -1.0);
		
		throw new IllegalArgumentException ("Unknown function " + fx);
	}
	
	/**
	 * This method returns all the operators defined in Polynomial.
	 * @return
	 */
	private static char[] getOperators () {
		return new char [] {'-', '+', '*', '/', '^', '%'};
	}
	
	/**
	 * This method returns all the operators defined in Polynomial.
	 * @return
	 */
	private static String[] getOperatorsAsString () {
		return new String [] {"-", "+", "*", "/", "^", "%"};
	}
	
	// method to check if a character is an operator
	private static boolean isOperator (char c)
	{
		for (char op : getOperators())
			if (c == op)
				return true;
		return false;
	}
	
	// method to check if a String is an operator
	private static boolean isOperator (String c) {
		for (char op : getOperators())
			if (c.equals("" + op))
				return true;
		return false;
	}
	
	// method to get the available functions
	private static String [] getFunctions()
	{
		return new String [] {"ln", "exp", "sin", "cos", "tan",
				"sinh", "cosh", "tanh", "cosec", "sec", "cot",
				"cosech", "sech", "coth"};
	}
	
	// method to check whether a function is available
	private static boolean isFunction (String fx)
	{
		for (String f : getFunctions())
			if (f.equals(fx))
				return true;
		return false;
	}
	
	// method that assigns a function node by mapping name to the fuction
	private static ExpNode assignFunctionNode (String fx, ExpNode operand)
	{	
		if (fx.equals("ln")) return new LnNode (operand);
		if (fx.equals("exp")) return new ExponentialNode (operand);
		if (fx.equals("sin")) return new SinNode(operand);
		if (fx.equals("cos")) return new CosNode (operand);
		if (fx.equals("tan")) return new TanNode (operand);
		if (fx.equals("cosec"))	return new CosecNode(operand);
		if (fx.equals("sec")) return new SecNode (operand);
		if (fx.equals("cot")) return new CotNode (operand);
		if (fx.equals("sinh")) return new SinhNode(operand);
		if (fx.equals("cosh")) return new CoshNode(operand);
		if (fx.equals("tanh")) return new TanhNode(operand);
		if (fx.equals("cosech")) return new CosechNode(operand);
		if (fx.equals("sech")) return new SechNode (operand);
		if (fx.equals("coth")) return new CothNode (operand);
		
		throw new IllegalArgumentException ("Unknown function " + fx);
	}
	
	// This method compare the precendence of the two operators.
	// True is returned when stackOp is equal or higher than infixOp or 
	// when stackOp is accidentally a function
	private static boolean precedence (char infixOp, String stackOp) 
			throws UnknownOperatorException
	{
		// check that op1 and op2 are both operators
		if (!(isOperator (infixOp) || isFunction ("" + infixOp)) &&
				!(isOperator (stackOp) || isFunction (stackOp)))
			throw new UnknownOperatorException (String.format(
					"Either operator '%s' or '%s' is invalid",
					infixOp, stackOp));
		
		// check for the case of equality
		if (stackOp.equals("" + infixOp))
			return true;
		
		// check for comparison with a function
		if (isFunction(stackOp))
			return true;
		
		// get the indices of the operators and compare
		int infixOpIndex = getOperatorHierachy ("" + infixOp);
		int stackOpIndex = getOperatorHierachy (stackOp);
		
		return stackOpIndex >= infixOpIndex;
	}
	
	// this method assigns an integer to operators in order to predict
	// their precedence.
	private static int getOperatorHierachy (String op)
	{
		if (op.equals("-") || op.equals("+"))
			return 1;
		else if (op.equals("*") || op.equals("/"))
			return 2;
		else if (op.equals("%") || op.equals("^"))
			return 3;
		else if (isFunction(op))
			return 4;
		throw new UnknownOperatorException ("Invalid operator " + op);
	}
	
	/**
	 * This method returns the root of a polynomial whose coefficients are the 
	 * element of the array a. The polynomial is of the form where the inital guess
	 * rr = 1, ss = 1, maxit = 10000, error = 0.000001
	 * 		f(x) = a0 + a1 x + a2 x^2 + ... + an-1 x^(n - 1) + an x^n = 0
	 * @see getRoots (a, n, rr, ss, maxit, error)
	 * @param a The coefficients of the polynomial
	 * @param n The order of the polynomial equation
	 * @return The roots as an array of complex numbers
	 * @throws IllegalArguementException when the length of array a does not correspond
	 * to the order of the Polynomial or the order is 0
	 */
	public static ComplexNumber[] getRoots (double [] a, int n)
	{
		return getRoots (a, n, 1, 1, 10000, 0.000001);
	}
	
	/**
	 * This method returns all the roots of a polynomial equation defined as
	 * 		f(x) = a0 + a1 x + a2 x^2 + ... + an-1 x^(n - 1) + an x^n = 0
	 * where the coefficients are the element of array a
	 * @param a The coefficients of the polynomial
	 * @param n The order of the polynomial equation
	 * @param rr The initial guess for r
	 * @param ss The initial guess for s
	 * @param maxit The maximum iteration
	 * @param error The minimum percentage error in the evaluation of r and s
	 * @return The roots as an array of complex numbers
	 * @throws IllegalArguementException when the length of array a does not correspond
	 * to the order of the Polynomial or the order is 0
	 */
	public static ComplexNumber[] getRoots (double [] a, int n, double rr, double ss,
			int maxit, double error)
	{
		// check for valid input
		if (a.length != n + 1)
			throw new IllegalArgumentException ("The length of a does not correspond " +
					"with the order of the polynomial");
		if (n == 0)
			throw new IllegalArgumentException ("The order of a polynomial cannot be" +
					" zero");
		
		// local variables
		double [] b = new double [n + 1];
		double [] c = new double [n + 1];
		ComplexNumber [] roots = new ComplexNumber[n];
		double r = rr, s = ss;
		double er1 = 1.0, er2 = 1.0;
		
		// iterate to evaluate the roots greater than 2
		while (n > 2) {
			// find the best values for r and s
			for (int i = 0; i < maxit; ++i) {
				b[n] = a[n];
				b[n - 1] = a[n - 1] + r * b[n];
				c[n] = b[n];
				c[n - 1] = b[n - 1] + r * c[n];
				
				for (int j = n - 2; j >= 0; --j) {
					b[j] = a[j] + r * b[j + 1] + s * b[j + 2];
					c[j] = b[j] + r * c[j + 1] + s * c[j + 2];
				}
				
				// evaluate the determinant for the alteration values of r and s
				double det = c[2] * c[2] - c[3] * c[1];
				if (det != 0) {
					double dr = (-b[1] * c[2] + b[0] * c[3]) / det;
					double ds = (-b[0] * c[2] + b[1] * c[1]) / det;
					r = r + dr;
					s = s + ds;
					
					// check the error
					if (r != 0)
						er1 = Math.abs(dr / r) * 100;
					if (s != 0)
						er2 = Math.abs(ds / s) * 100;
				} else {
					r = r + 1;		// readjust the initial guess value for r and s
					s = s + 1;
					i = 0;			// restart the iteration for values of r and s
				}
				
				// check the error
				if (er1 <= error && er2 <= error)
					break;
			}	// end the loop for evaluating values for r and s
			
			// get the quadratic roots
			ComplexNumber[] newRoots = getQuadraticRoot (r, s);
			for (int i = 1; i < 3; ++i)
				roots[n - i] = newRoots[i - 1];
			
			// revaluate the values of n and a for the next iterate
			n -= 2;
			for (int j = 0; j <= n; ++j)
				a[j] = b[j + 2];
		}	// end the loop for evaluation of roots
		
		// evaluate the rest
		if (n == 2) {
			r = -a[1] / a[2];
			s = -a[0] / a[2];
			
			// get the quadratic roots
			ComplexNumber[] newRoots = getQuadraticRoot (r, s);
			for (int i = 1; i >= 0; --i)
				roots[i] = newRoots[i];
		} else if (n == 1)
			roots[0] = new ComplexNumber(-a[0] / a[1], 0);
		
		return roots;
	}
	
	/**
	 * This method evaluate the root of a quadratic equation of the for
	 * 		ax^2 - rx - s = 0
	 * where a = 1
	 * @param r
	 * @param s
	 * @return
	 */
	private static ComplexNumber [] getQuadraticRoot (double r, double s)
	{
		// local variable
		ComplexNumber [] root = new ComplexNumber [2];
		
		// evaluate the discreminant
		double disc = r * r + 4 * s;
		
		if (disc > 0) {
			root[0] = new ComplexNumber ((r + Math.sqrt(disc)) / 2.0);
			root[1] = new ComplexNumber ((r - Math.sqrt(disc)) / 2.0);
		} else {
			root[0] = new ComplexNumber (r / 2.0, (Math.sqrt(-disc)) / 2.0);
			root[1] = root[0].getConjugate();
		}
		return root;
	}
}