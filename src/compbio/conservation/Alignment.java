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
}
