package compbio.conservation;

public class ColumnCollection {

	private final Column[] cols;
	
	private double[] kabat = null;
	
	private double[] jores = null;

	private double[] schneider = null;

	private double[] shenkin = null;

	private double[] gerstein = null;

	private double[] taylorNoGaps = null;

	private double[] taylorGaps = null;

	private double[] zvelibil = null;

	private double[] karlin = null;

	private double[] armon = null;

	private double[] thompson = null;

	private double[] lancet = null;

	private double[] mirny = null;

	private double[] williamson = null;

	private double[] landgraf = null;

	private double[] sander = null;

	private double[] valdar = null;

	public ColumnCollection(AminoAcidMatrix m) {
	
	if (m == null) {
		
		throw new IllegalArgumentException("m must not be null");
	}
	
	assert m.numberOfColumns() != 0 : "Sth wrong with the matrix, not null but doesn't have any residues";
	
	cols = new Column[m.numberOfColumns()];
	
	for (int i = 0; i < m.numberOfColumns(); i++) {
		
		cols[i] = new Column(m, i);
		
	}
	
	}	
	
	Column[] getColumnCollection () {
	
	return cols;
	
	}

	int collectionLength() {
	
	int len = cols.length;
	
	return len;
	
	}

	void calculateKabat() {
		
		kabat = new double[cols.length];
		
		for (int i = 0; i < cols.length; i++) {
			
			kabat[i] = cols[i].kabatScore();
		}
	
	}
	
	double[] getKabat() {
		
		if ( kabat == null) {
			
			this.calculateKabat();
		}
		
		return kabat;
		
	}
	
	void calculateSchneider() {
		
		schneider = new double[cols.length];
		
		for (int i = 0; i < cols.length; i++) {
			
			schneider[i] = cols[i].schneiderScore();
		}
	
	}
	
	double[] getSchneider() {
		
		if ( schneider == null) {
			
			this.calculateSchneider();
		}
		
		return schneider;
		
	}
	
	void calculateShenkin() {
		
		shenkin = new double[cols.length];
		
		for (int i = 0; i < cols.length; i++) {
			
			shenkin[i] = cols[i].schneiderScore();
		}
	
	}
	
	double[] getShenkin() {
		
		if ( shenkin == null) {
			
			this.calculateShenkin();
		}
		
		return shenkin;
		
	}

	void calculateGerstein() {
		
		gerstein = new double[cols.length];
		
		for (int i = 0; i < cols.length; i++) {
			
			gerstein[i] = cols[i].gersteinScore();
		}
	}
	
	double[] getGerstein() {
		
		if ( gerstein == null) {
			
			this.calculateSchneider();
		}
		
		return schneider;
		
	}

	void calculateTaylorNoGaps() {
		
		taylorNoGaps = new double[cols.length];
		
		for (int i = 0; i < cols.length; i++) {
			
			taylorNoGaps[i] = cols[i].SmallestTaylorSetNoGaps();
		}
	
	}
	
	double[] getTaylorNoGaps() {
		
		if ( taylorNoGaps == null) {
			
			this.calculateTaylorNoGaps();
		}
		
		return taylorNoGaps;
		
	}
	
	void calculateTaylorGaps() {
		
		taylorGaps = new double[cols.length];
		
		for (int i = 0; i < cols.length; i++) {
			
			taylorGaps[i] = cols[i].SmallestTaylorSetGaps();
		}
	
	}
	
	double[] getTaylorGaps() {
		
		if ( taylorGaps == null) {
			
			this.calculateTaylorGaps();
		}
		
		return taylorGaps;
		
	}

	void calculateZvelibil() {
		
		zvelibil = new double[cols.length];
		
		for (int i = 0; i < cols.length; i++) {
			
		//	zvelibil[i] = cols[i].zvelibilScore();
		}
	
	}
	
	double[] getZvelibil() {
		
		if ( zvelibil == null) {
			
			this.calculateZvelibil();
		}
		
		return zvelibil;
		
	}

	void calculateKarlin() {
		
		karlin = new double[cols.length];
		
		for (int i = 0; i < cols.length; i++) {
			
			karlin[i] = cols[i].karlinScore();
		}
	
	}
	
	double[] getKarlin() {
		
		if ( karlin == null) {
			
			this.calculateKarlin();
		}
		
		return karlin;
		
	}

	void calculateArmon() {
		
		armon = new double[cols.length];
		
		for (int i = 0; i < cols.length; i++) {
			
			armon[i] = cols[i].armonScore();
		}
	
	}
	
	double[] getArmon() {
		
		if ( armon == null) {
			
			this.calculateArmon();
		}
		
		return armon;
		
	}

	void calculateThompson() {
		
		thompson = new double[cols.length];
		
		for (int i = 0; i < cols.length; i++) {
			
			thompson[i] = cols[i].thompsonScore();
		}
	
	}
	
	double[] getThompson() {
		
		if ( thompson == null) {
			
			this.calculateThompson();
		}
		
		return thompson;
		
	}

	void calculateLancet() {
		
		lancet = new double[cols.length];
		
		for (int i = 0; i < cols.length; i++) {
			
			lancet[i] = cols[i].lancetScore();
		}
	
	}
	
	double[] getLancet() {
		
		if ( lancet == null) {
			
			this.calculateLancet();
		}
		
		return lancet;
		
	}

	void calculateMirny() {
		
		mirny = new double[cols.length];
		
		for (int i = 0; i < cols.length; i++) {
			
			mirny[i] = cols[i].mirnyScore();
		}
	
	}
	
	double[] getMirny() {
		
		if ( taylorGaps == null) {
			
			this.calculateMirny();
		}
		
		return mirny;
		
	}

	void calculateWilliamson() {
		
		williamson = new double[cols.length];
		
		for (int i = 0; i < cols.length; i++) {
			
			williamson[i] = cols[i].williamsonScore();
		}
	
	}
	
	double[] getWilliamson() {
		
		if ( williamson == null) {
			
			this.calculateWilliamson();
		}
		
		return williamson;
		
	}

	void calculateLandgarf() {
		
		landgraf = new double[cols.length];
		
		for (int i = 0; i < cols.length; i++) {
			
			landgraf[i] = cols[i].landgrafScore();
		}
	
	}
	
	double[] getLandgraf() {
		
		if ( landgraf == null) {
			
			this.calculateLandgarf();
		}
		
		return landgraf;
		
	}

	void calculateSander() {
		
		sander = new double[cols.length];
		
		for (int i = 0; i < cols.length; i++) {
			
			sander[i] = cols[i].sanderScore();;
		}
	
	}
	
	double[] getSander() {
		
		if ( sander == null) {
			
			this.calculateSander();
		}
		
		return sander;
		
	}
	
	void calculateValdar() {
		
		valdar = new double[cols.length];
		
		for (int i = 0; i < cols.length; i++) {
			
			valdar[i] = cols[i].valdarScore();
		}
	
	}
	
	double[] getValdar() {
		
		if ( valdar == null) {
			
			this.calculateValdar();
		}
		
		return valdar;
		
	}

	void calculateJores() {
		
		jores = new double[cols.length];
		
		for (int i = 0; i < cols.length; i++) {
			
			jores[i] = cols[i].joresScore();
			
		}
		
	}
	
	double[] getJores() {
		
		if( jores == null) {
			
			this.calculateJores();
		}
	
		return jores;
	}
	
}





