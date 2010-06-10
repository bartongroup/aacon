package compbio.conservation;

import java.util.*;

public class ColumnCollection {

private Column[] col;

public ColumnCollection(AminoAcidMatrix m) {
	
	if (m == null) {
		
		throw new IllegalArgumentException("m must not be null");
	}
	
	assert m.numberOfColumns() != 0 : "Sth wrong with the matrix, not null but doesn't have any residues";
	
	col = new Column[m.numberOfColumns()];
	
	for (int i = 0; i < m.numberOfColumns(); i++) {
		
		col[i] = new Column(m, i);
		
	}
	
}	
	
Column[] getColumnCollection () {
	
	return col;
}

int collectionLength() {
	
	int len = col.length;
	
	return len;
}
// this method calculates the overall amino acid frequency over the column window present
// still to be written

Map<Character,Integer> TotalAcids() {
	
	Map<Character,Integer>  freq = new HashMap<Character,Integer>();
	
	
	return freq;
	
}
}



