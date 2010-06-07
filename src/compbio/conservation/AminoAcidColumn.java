package compbio.conservation;

//this is class Column, the objects created represent columns of an alignment
//column collection will be created
//columns are than used in various calculations
//the idea is that all the columns objects will be stored in an array
//that array will be created as instantiation of ColumnCollection class
//will see what's gonna happen

public class AminoAcidColumn {
	
	private char[] columnArr;
	private AminoAcidOccurance[] aaOcc;

	// constructor, 

	public AminoAcidColumn( AminoAcidMatrix m, int c) {

	int colLength = m.numberOfRows();

	StringBuilder s = new StringBuilder();

	columnArr = new char[colLength];

	// creates a string 

	for ( int i = 0; i < colLength ; i++) {
		
	        char aacid = m.getMatrixPosition( i, c );
		s = s.append(aacid);
	}
	String columnStr = s.toString().toUpperCase();
	
	columnArr = columnStr.toCharArray();

	// creates occurance objects
	
	this.createOccurance();
	
	// counts occurance
	
	this.countOccurance();
	 
	}

    //returns column length

	public int length(){

	int len = columnArr.length;

	return len;

	}

	// checks whether column is empty(consists of gaps only)

	public boolean isEmpty() {

	int nrOfGaps = 0;

	int len = columnArr.length;

	for (int i = 0; i < len; i++) {

	if (columnArr[i] == '-')
	      nrOfGaps++;
	}

	if (nrOfGaps == len)
	      return true;

	else
	      return false;

	}

	// checks whether there are gaps present in the column,
	//returns true if there are, false if there aren't any

	public boolean gapsInColumn() {

	int nrOfGaps = 0;

	int len = columnArr.length;

	for (int i = 0; i < len; i++) {

	if (columnArr[i] == '-')
	      nrOfGaps++;
	}

	if (nrOfGaps > 0)
	      return true;

	else
	      return false;

	}

	// retuns true if all but one residues in the column are gaps,
	// false otherwise

	public boolean allButOneGaps() {

	int nrOfGaps = 0;

	int len = columnArr.length;

	for (int i = 0; i < len; i++) {

	if (columnArr[i] == '-')
	      nrOfGaps++;
	}

	if (nrOfGaps == len - 1)
	      return true;

	else
	      return false;

	}

	// count the numer how many times each aa occurs in the in the column
	// the information is stored in an array 
	// of AminoAcidOccuranceChars objects that were already instantiated
	// thinking about making it a private method
	// and moving to the constructor

	private void countOccurance() { 
		
	AminoAcidAlphabet alp = new AminoAcidAlphabet(); 
	
	char aAcids[] = alp.getAlphabet();
	
    for( int i = 0; i < columnArr.length; i++) {

	     for( int j = 0; j < aAcids.length; j++) {

	           if (columnArr[i] == aaOcc[j].getId()) {

	                     aaOcc[j].addToOccurance(1);

	                     break;
	              }

	        }

	   }

	}

	// calculates how many residue types are  present
    public int howManyResidueTypes() { 
		
    int resTypeNr = 0;	

	for ( int i = 0; i < aaOcc.length; i++ ) {
		
		if(aaOcc[i].getOccurance() > 0) 
			resTypeNr++;
	}
	
		return resTypeNr;
		
	}
    
  // creates occurance objects, sets occurance to 0
  // creates 21 AminoAcidOccurance objects, the occurance is set to 0
  // id is set to one of 20 aa(one letter code) or gap(dash)

    
    private void createOccurance() {
    	
    	AminoAcidAlphabet alp = new AminoAcidAlphabet();
    	
    	char[] aAcids = alp.getAlphabet();
    	
    	aaOcc = new AminoAcidOccurance[alp.length()];

    	for ( int i = 0; i < alp.length(); i++ ) {

    	aaOcc[i] = new AminoAcidOccurance( aAcids[i]);

    	}
    	
    	}
    
    // gives the most common aa occurance 
    
    public int  mostCommonOccurance() {
    	
    int mco = 0;
    
    for (int i = 0; i < aaOcc.length; i++ ) {
    	
    	if (aaOcc[i].getOccurance() > mco)
    		 mco = aaOcc[i].getOccurance();
    	
    }
    
    return mco;
    
    }
    	
    }
    

	

