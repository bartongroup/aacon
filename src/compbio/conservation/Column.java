package compbio.conservation;
import java.util.Map;
import java.util.HashMap;
import org.testng.Assert;
import java.util.Collection;
import java.util.Iterator;




public class Column {

	private final char[] columnArr;
    private final Map<Character,Integer> acidsIntMap;

    // TO BE DELETED, THIS CONSTRUCTOR USED FOR TESTS ONLY

    //public Column(char a, char b, char c, char d, char e, char f,
	 //  char g, char h, char i, char j) {

	//columnArr = new char[10];

	//columnArr[0] = a;
	//columnArr[1] = b;
	//columnArr[2] = c;
	//columnArr[3] = d;
	//columnArr[4] = e;
	//columnArr[5] = f;
	//columnArr[6] = g;
	//columnArr[7] = h;
	//columnArr[8] = i;
	//columnArr[9] = j;


    //}

    // constructor, 

    public Column(AminoAcidMatrix m, int column) {

	if (m == null) {
	    throw new IllegalArgumentException("Matrix must not be null");
	}
	assert column >= 0 && column <= m.numberOfColumns();

	int colLength = m.numberOfRows();

	StringBuilder s = new StringBuilder();

	// creates a string 

	for (int i = 0; i < colLength; i++) {

	    char aacid = m.getMatrixPosition(i, column);
	    s = s.append(aacid);
	}

	String columnStr = s.toString().toUpperCase();

	columnArr = columnStr.toCharArray();

// calls mathod to create a map
	
	Alphabet alp = new Alphabet();
	
	acidsIntMap = alp.calculateOccurance(columnArr);
	
	}
	
    // check whether column is empty 
    
    boolean isEmpty() {
    	
    	if (acidsIntMap.containsKey('-') && acidsIntMap.get('-') == columnArr.length) {
    			
    			return true;
    		}
    	
    	else {
    		
    		return false;
    		
    	}
    	
    }
	// returns true if all but one residues in the column are gaps
    boolean allButOneGaps() {
    	
    	if (acidsIntMap.containsKey('-') && acidsIntMap.get('-') == columnArr.length - 1) {
			
			return true;
		}
	
	else {
		
		return false;
		
	}
	
    }
    // returns true if there is only one residue type in the column
    boolean oneResidueType() {
    	
    	if (acidsIntMap.size() == 1 && acidsIntMap.containsKey('-') == false) {
    		
    		return true;
    		
    	}
    	
    	else {
    		
    		return false;
    		
    	}
    	
    }
 
    // returns the nr of aa types in the column
    int numberOfAcids() {
    	
    	return acidsIntMap.size(); 
    	
    }
    
    int mostCommonNumber() {
    	
    	int max = 0;
    	
    	Collection<Integer> values = acidsIntMap.values();
    	
    	Iterator<Integer> itr = values.iterator();
    	
    	while(itr.hasNext()) {
    		
    		int val = itr.next();
    		
    		if (val > max) {
    			
    			max = val;
    			
    		}
    		
    	}
    	
    	assert max != 0 : "Zero in the most Common Number";
    	
    	return max;
    	
    }
    
    Map<Character,Integer> getAcidsIntMap() {
    
    return acidsIntMap;
    
    
    }
    			
 int length() {
	 
	 return columnArr.length;
	 
 }
 }
  
    	
    	
 

