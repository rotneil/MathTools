package phDModels;

import math.tools.ExpNode;
import math.tools.Matrix;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Nehemiah Oluwafemi
 */
public class ModelDialog extends javax.swing.JDialog {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * The dialog listener.
	 * @author Nehemiah Oluwafemi
	 */
	public interface ModelDialogListener {
		/**
		 * The method that is called when a new model has just been initialized.
		 * 
		 * @param modelName This is the name that the system is called
		 * @param fx The State equation as an array of ExpNode
		 * @param x The independent coordinates of the model
		 * @param variables The state variables including the independent 
		 * coordinates x and the parameters
		 * @param values The corresponding initial values for the variables
		 */
		public void onModelInitialized (String modelName, ExpNode[] fx, String [] x, 
				String [] variables, double [] values, int [] coord);
	}
	
	// instance variables
	private ModelDialogListener mListener;
	javax.swing.JFrame mFrame;
	
    /**
     * This creates the dialog that allows the user to enter a state model
     * @param mFrame The calling frame
     * @param l The ModelDialogListener
     */
    public ModelDialog (javax.swing.JFrame frame, ModelDialogListener l) {
        super(frame, "State Model Initializer Dialog", true);
        initComponents();
        
        mListener = l;
        mFrame = frame;
        setLocationRelativeTo (frame);
    }                    

    private void cancleButtonActionPerformed(java.awt.event.ActionEvent evt) {                                             
        // TODO add your handling code here:
    	// check that the user has not entered anything
    	if (modelArea.getText().isEmpty())
    		setVisible (false);
    	else {
    		int response = javax.swing.JOptionPane.showConfirmDialog(
    				this, "Are sure you want to cancel? ", 
    				"POSSIBLE LOSS OF DATA",
    				javax.swing.JOptionPane.YES_NO_OPTION);
    		if (response == javax.swing.JOptionPane.YES_OPTION)
    			setVisible (false);
    	}
    }                                            

    private void enterButtonActionPerformed(java.awt.event.ActionEvent evt) {                                            
        // TODO add your handling code here:
    	try {
    		// process user's input
    		String modelName = nameField.getText();
    		if (modelName.isEmpty())
    			throw new IllegalArgumentException ("The name of the model " +
    					"cannot be empty");
        	String modelInput = modelArea.getText();
        	String varInput = variableField.getText();
        	String valInput = valueField.getText();
        	
        	// convert the entered variables and values into array
        	String [] var = varInput.split(",");
        	String [] valStr = valInput.split(",");
        	String [] coordStr = coordField.getText().split(",");
        	
        	// confirm equal length of variables and values
        	if (var.length != valStr.length)
        		throw new IllegalArgumentException ("The number of Variables " +
        				"does not correspond to the number of Values");
        	
        	// throw exception if the plot coordinate is more than 2
        	if (coordStr.length != 2)
        		throw new IllegalArgumentException ("The plot coordinate " +
        				"can only be two");
        	
        	// trim white spaces off var
        	for (int i = 0; i < var.length; ++i) {
        		while (var[i].startsWith(" "))
        			var[i] = var[i].substring(1, var[i].length());
        		while (var[i].endsWith(" "))
        			var[i] = var[i].substring(0, var[i].length() - 1);
        	}
        	
        	// convert values to double
        	double [] val = new double [var.length];
        	for (int i = 0; i < var.length; ++i)
        		val[i] = Double.parseDouble(valStr[i]);
        	
        	// extract the model equation from each modelInput
        	String [] model = modelInput.split("\n");
        	java.util.ArrayList<StringBuffer> eq = 
        			new java.util.ArrayList<StringBuffer>();
        	java.util.ArrayList<String> x = new java.util.ArrayList<String>();
        	String [] op;
        	for (int i = 0; i < model.length; ++i) {
        		// leave model if it is just an empty line
        		if (!model[i].contains("="))
        			continue;
        		
        		// split arround
        		op = model[i].split("=");
        		
        		// check if left hand has something like 'd(x1)/dt'
        		if (op[0].contains("("))
        			x.add(op[0].substring (op[0].indexOf("(") + 1, 
        					op[0].indexOf(")")));
        		else
        			x.add(op[0]);
        		
        		// initialize equation
        		eq.add(new StringBuffer (op[1]));
        	}
        	StringBuffer [] eqBuffer = new StringBuffer [eq.size()];
        	String [] xx = new String [x.size()];
        	eq.toArray(eqBuffer);
        	x.toArray(xx);
        	
        	// confirm that the equation is consistent
        	if (eqBuffer.length != xx.length)
        		throw new IllegalArgumentException ("Inconsistency in state " +
        				"equation definition");
        	        	
        	// trim off whitespaces from coord and then get the indices
        	int [] coord = new int [coordStr.length];
        	for (int i = 0; i < coord.length; ++i) {
        		while (coordStr[i].startsWith(" "))
        			coordStr[i] = coordStr[i].substring(1, coordStr[i].length());
        		while (coordStr[i].endsWith(" "))
        			coordStr[i] = coordStr[i].substring(0, coordStr[i].length() - 1);
        		coord[i] = ExpNode.getIndex(coordStr[i], xx);
        	}
        	
        	// convert eq to ExpNode
        	ExpNode [] fx = Matrix.infixToExpNode(eqBuffer, var);
        	
        	// prompt listener for action
        	if (mListener != null)
        		mListener.onModelInitialized(modelName, fx, xx, var, val, coord);
        	
        	// show success message and dispose dialog
        	setVisible (false);
        	javax.swing.JOptionPane.showMessageDialog(mFrame,
        			"Model initialized successfully", "SUCCESS MESSAGE",
        			javax.swing.JOptionPane.INFORMATION_MESSAGE);
    	} catch (IllegalArgumentException e) {
    		javax.swing.JOptionPane.showMessageDialog(this, 
    				e.getMessage(), e.getClass().getSimpleName(),
    				javax.swing.JOptionPane.ERROR_MESSAGE);
    	} catch (NullPointerException e) {
    		javax.swing.JOptionPane.showMessageDialog(this, 
    				e.getMessage(), e.getClass().getSimpleName(),
    				javax.swing.JOptionPane.ERROR_MESSAGE);
    	}
    }

    // Variables declaration - do not modify                     
    private javax.swing.JButton cancelButton;
    private javax.swing.JButton enterButton;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JTextArea modelArea;
    private javax.swing.JTextField valueField;
    private javax.swing.JTextField variableField;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JTextField coordField;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JTextField nameField;
    
    // End of variables declaration
    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    // <editor-fold defaultstate="collapsed" desc="Generated Code">                          
    private void initComponents() {
        jPanel1 = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        jLabel2 = new javax.swing.JLabel();
        variableField = new javax.swing.JTextField();
        jScrollPane1 = new javax.swing.JScrollPane();
        modelArea = new javax.swing.JTextArea();
        jLabel3 = new javax.swing.JLabel();
        valueField = new javax.swing.JTextField();
        jLabel4 = new javax.swing.JLabel();
        coordField = new javax.swing.JTextField();
        cancelButton = new javax.swing.JButton();
        enterButton = new javax.swing.JButton();
        jLabel5 = new javax.swing.JLabel();
        nameField = new javax.swing.JTextField();
        
        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);

        jLabel1.setText("Enter the State Equation");

        jLabel2.setText("Enter the State Variables as a comma seperated strings");

        modelArea.setColumns(20);
        modelArea.setLineWrap(true);
        modelArea.setRows(5);
        modelArea.setWrapStyleWord(true);
        jScrollPane1.setViewportView(modelArea);

        jLabel3.setText("Enter the corresponding initial Values for the state variables");
        
        jLabel4.setText("Enter the plot coordinates x,y as comma separated variables");
        
        jLabel5.setText("Enter the model name for the system");
        
        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cancleButtonActionPerformed(evt);
            }
        });

        enterButton.setText("Enter");
        enterButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                enterButtonActionPerformed(evt);
            }
        });

        coordField.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                enterButtonActionPerformed(evt);
            }
        });
        
        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
            		.addComponent(jLabel5, javax.swing.GroupLayout.DEFAULT_SIZE, 305, Short.MAX_VALUE)
                    .addComponent(nameField)
                    .addComponent(jLabel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jScrollPane1)
                    .addComponent(jLabel2, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(variableField)
                    .addComponent(jLabel3, javax.swing.GroupLayout.DEFAULT_SIZE, 305, Short.MAX_VALUE)
                    .addComponent(valueField)
                    .addComponent(jLabel4, javax.swing.GroupLayout.DEFAULT_SIZE, 305, Short.MAX_VALUE)
                    .addComponent(coordField)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(cancelButton, javax.swing.GroupLayout.PREFERRED_SIZE, 151, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(2, 2, 2)
                        .addComponent(enterButton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel5)
                .addGap(2, 2, 2)
                .addComponent(nameField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(2, 2, 2)
                .addComponent(jLabel1)
                .addGap(2, 2, 2)
                .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jLabel2)
                .addGap(2, 2, 2)
                .addComponent(variableField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(2, 2, 2)
                .addComponent(jLabel3)
                .addGap(2, 2, 2)
                .addComponent(valueField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addComponent(jLabel4)
                .addGap(2, 2, 2)
                .addComponent(coordField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(cancelButton)
                    .addComponent(enterButton))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
        );
        
        pack();
        
        setResizable (false);
    }
}