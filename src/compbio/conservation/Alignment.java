package compbio.conservation;

import java.util.*;

public class Alignment {
	
	private ColumnCollection cols;
	
	
	public Alignment (AminoAcidMatrix m){
		
		if (m == null) {
			
			throw new IllegalArgumentException("matrix must not be zero");
		}
		
		
	cols = new ColumnCollection(m);
	
	
	}

	
double[] kabat() {
	 
	double[] result = new double[cols.getColumnCollection().length];
	
	Column[] columns = cols.getColumnCollection();
	
	assert columns.length > 0 ;
	
	for (int i = 0; i < result.length; i++ ) {
		
		if (columns[i].isEmpty() == true) {
			
			result[i] = 0;
			
		}
	
		else {
			
		result[i] = columns[i].length() * columns[i].numberOfAcidsNoGap()/ columns[i].mostCommonNumber(); 
		
		}
	}
	
	return result;
	
}

// Symbol Enthropy Scores

double[] schneider() {
	
	double[] result = new double[cols.getColumnCollection().length]; // might change for matrix length if matrix will get a set size
	
	double normal = 1.0 / 20.0;
	
	Column[] columns = cols.getColumnCollection();
	
	assert columns.length > 0 ;
	
	for (int i = 0; i < result.length; i++) {
		
		result[i] = ShannonEnthropy.ShannonLn(columns[i].getAcidsIntMap(), columns[i].length()) * normal;
		
	}

	return result;
}


double[] shenkin() {
	
	double[] result = new double[cols.getColumnCollection().length]; 
	
	Column[] columns = cols.getColumnCollection();
	
	assert columns.length > 0 ;
	
	for (int i = 0; i < result.length; i++) {
		
		result[i] = Math.pow( 2.0, ShannonEnthropy.ShannonLog2(columns[i].getAcidsIntMap(), columns[i].length())) * 6.0;
		
	}

	return result;
}

double[] gerstein() {
	
	double[] result = new double[cols.getColumnCollection().length]; 
	
	Column[] columns = cols.getColumnCollection();
	
	assert columns.length > 0 ;
	
	for (int i = 0; i < result.length; i++) {
		
		result[i] = ShannonEnthropy.ShannonLn(cols.TotalAcids(), cols.getColumnCollection().length * columns[0].length() ) - ShannonEnthropy.ShannonLn(columns[i].getAcidsIntMap(), columns[i].length());
		 
	}
	
	return result;
	
}

// stereochemical poperty scores


// have to figure out where to get the sets from
// initially create an object with sets, maybe something else, but the important part 
// is that we need a hashmap of sets
// sets used here might actually get stored as the instance variable of alignmanet class
// will see, ahve to think about it
// thoght about it, will create a text file with sets and a class and a staic method
// that reads in the sets an dputs thatm in hashmaps

int[] taylor() {
	
int[] result = new int[cols.getColumnCollection().length]; 
	
	Column[] columns = cols.getColumnCollection();
	
	assert columns.length > 0 ;
	
	Map<String, HashSet<Character>> set = new HashMap<String, HashSet<Character>>();
	
	for (int i = 0; i < result.length; i++) {
		
		//result[i] = columns[i].SmallestTaylorSet(set);
		
	}
	
	return result;
	
	}

}

