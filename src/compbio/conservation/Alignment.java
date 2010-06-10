package compbio.conservation;

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
		
		result[i] = columns[i].length() * columns[i].numberOfAcids()/ columns[i].mostCommonNumber(); 
		
	}
	
	return result;
	
}

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





}

