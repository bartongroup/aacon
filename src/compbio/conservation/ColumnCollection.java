package compbio.conservation;

public class ColumnCollection {

private AminoAcidColumn[] col;

public ColumnCollection(AminoAcidMatrix m) {
	
	col = new AminoAcidColumn[m.numberOfColumns()];
	
	for (int i = 0; i < m.numberOfColumns(); i++) {
		
		col[i] = new AminoAcidColumn(m, i);
		
	}
	
}	
	
public AminoAcidColumn[] getColumnCollection () {
	
	return col;
}

public int collectionLength() {
	
	int len = col.length;
	
	return len;
}

}
