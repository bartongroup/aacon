package compbio.conservation;

public class Alignment {
	
	private AminoAcidColumn[] cols;
	
	public Alignment (AminoAcidMatrix m){
		
	ColumnCollection col = new ColumnCollection(m);
	
	cols = col.getColumnCollection();
	
	}

	
public double[] kabat() {
	 
	double[] result = new double[cols[0].length()];
	
	for (int i = 0; i < result.length; i++ ) {
		
		result[i] = cols[i].length() * cols[i].howManyResidueTypes() / cols[i].length(); 
		
	}
	
	return result;
	
}

public double[] jores() {
	
	double[] result = new double[cols[0].length()];
	
	for (int i = 0; i < result.length; i++) {
		
		if (cols[i].isEmpty() == true)
			result[i] = 0;
		else {
		
		if (cols[i].allButOneGaps() == true)
			result[i] = 0;
		else {
		
		if (cols[i].onlyOneResType() == true)
			result[i] = 1;
		
		else { 
		
			PairCollection pair = cols[i].pairs2();
		
		    result[i] = (pair.getNrPairs() / pair.getMostFrequent()) * (cols[i].length() * (cols[i].length() - 1) / 2);
		    
		}
	
	    }
	
	    }

        }
	
	return result;
	
	}

}
