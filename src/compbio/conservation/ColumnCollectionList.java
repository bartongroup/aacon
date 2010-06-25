package compbio.conservation;

import java.util.*;

/**
 * Class that provides a framework for a container of many column objects.
 * It behaves like a window over a small part of of the alignment.
 * In this class, there are no window overlaps.
 * The window is "jumps", first it covers the n first columns, than it is moved and covers the next n of columns.
 *  
 * @author agolicz
 *
 */

public class ColumnCollectionList {
		/**
		 * Keeps track of how many changes of columns took place.
		 */
		private int changesOfColumns = 0;
	
		/**
		 * Keeps track of the number of   updates
		 */
		
		private int nrOfUpdates = 0;
		/**
		 * Holds the columns currently present in the window.
		 */
		private final Column[] cols;
		
		/**
		 * Holds information whether calculations were initialized after constructor was called. 
		 */
		
		private boolean calcsInitialized = false;
		
		/**
		 * Holds the list of Kabat results.
		 */
		
		private List<Double> kabat = null;
		
		/**
		 * Kabat initialization status
		 */
		
		private boolean kabatInitStat = false;
		
		/**
		 * Holds list of Jores results.
		 */
		
		private List<Double> jores = null;
		
		/**
		 * Jores initialization status.
		 */
		
		private boolean joresInitStat = false;
		
		/**
		 * Holds list of Schneider results.
		 */

		private List<Double> schneider = null;
		
		/**
		 * Schneider initialization status.
		 */
		
		private boolean schneiderInitStat = false;
		
		/**
		 * Holds list of Shenkin results.
		 */

		private List<Double> shenkin = null;
		
		/**
		 * Shenkin initialization status.
		 */
		
		private boolean shenkinInitStat = false;
		
		/**
		 * Holds list of Gerstein results.
		 */

		private List<Double> gerstein = null;
		
		/**
		 * Gerstein initialization status.
		 */
		
		private boolean gersteinInitStat = false;
		
		/**
		 * Holds list of Taylor results.
		 */

		private List<Integer> taylorNoGaps = null;
		
		/**
		 * Taylor initialization status.
		 */
		
		private boolean taylorNoGapsInitSatat = false;
		
		/**
		 * Holds list of Taylor results.
		 */

		private List<Integer> taylorGaps = null;
		
		/**
		 * Taylor initialization status.
		 */
		
		private boolean taylorGapsInitStat = false;
		
		/**
		 * Holds list of Zvelibil results.
		 */

		private List<Integer> zvelibil = null;
		
		/**
		 * Zvelibil initialization status.
		 */
		
		private boolean zvelibilInitStat = false;
		
		/**
		 * Holds list of Karlin results.
		 */

		private List<Double> karlin = null;
		
		/** 
		 * Karlin initialization status.
		 */
		
		private boolean karlinInitStat = false;
		
		/**
		 * Holds list of Armon results.
		 */

		private List<Double> armon = null;
		
		/**
		 * Armon initialization status.
		 */
		
		private boolean armonInitStat = false;
		
		/**
		 * Holds list of Thompson Results.
		 */

		private List<Double> thompson = null;
		
		/**
		 * Thompson initialization status.
		 */
		
		private boolean thompsonInitStat = false;
		
		/**
		 * Holds list of Lancet results.
		 */

		private List<Double> lancet = null;
		
		/**
		 * Lancet initialization status.
		 */
		
		private boolean lancetInitStat = false;
		
		/**
		 * Holds list of Mirny results.
		 */

		private List<Double> mirny = null;
		
		/**
		 * Mirny initialization status.
		 */
		
		private boolean mirnyInitStat = false;
		
		/**
		 * List of Williamson results.
		 */

		private List<Double> williamson = null;
		
		/**
		 * Williamson initialization status.
		 */
		
		private boolean williamsonInitStat = false;
		
		/**
		 * List of landgraf results.
		 */

		private List<Double> landgraf = null;
		
		/**
		 * Landgraf initialization status.
		 */
		
		private boolean landgrafInitStat = false;
		
		/**
		 * List of Sander results.
		 */

		private List<Double> sander = null;
		
		/**
		 * Sander initialization status.
		 */
		
		private boolean sanderInitStat = false;
		
		/**
		 * List of Valdar Results.
		 */

		private List<Double> valdar = null;
		
		/**
		 * Valdar initialization Status.
		 */
		
		private boolean valdarInitStat = false;
		
		/**
		 * Holds reference to the "mother" matrix.
		 */
		
		private AminoAcidMatrix matrix;
		
		/**
		 * Constructor. Takes reference to matrix and window width as parameter.
		 * @param m
		 * @param nrOfColumns
		 */

		
		public ColumnCollectionList(AminoAcidMatrix m, int nrOfColumns) {
		
		if (m == null) {
			
			throw new IllegalArgumentException("m must not be null");
		}
		
		assert m.numberOfColumns() != 0 : "Sth wrong with the matrix, not null but doesn't have any residues";
		
		assert nrOfColumns < m.numberOfColumns();
		
		matrix = m;
		
		cols = new Column[nrOfColumns];
		
		for (int i = 0; i < nrOfColumns; i++) {
			
			cols[i] = new Column(matrix, i);
			
		}
		
		}	
		
		/**
		 * Returns columns currently present in the collection.
		 * 
		 * @return column collection
		 */
		
		Column[] getColumnCollection () {
		
		return cols;
		
		}
		
		/**
		 * Returns collection lenght( window width).
		 * 
		 * @return collection length
		 */

		int collectionLength() {
		
		int len = cols.length;
		
		return len;
		
		}
		
		/**
		 * Calculates Kabat for all columns in the collection.
		 */

		void calculateKabat() {
			
			for (int i = 0; i < cols.length; i++) {
				
				kabat.add(cols[i].kabatScore());
			}
		
		}
		
		/**
		 * Returns results for all the column there were and are present in the collection
		 * 
		 * @return results
		 */
		
		List<Double> getKabat() {
			
			
			assert this.nrOfUpdates == this.changesOfColumns && this.kabatInitStat == true;
			
			return kabat;
			
		}
		
		/**
		 * Calculates Schneider for all columns in the collection.
		 */
		
		void calculateSchneider() {
			
			
			for (int i = 0; i < cols.length; i++) {
				
				schneider.add(cols[i].schneiderScore());
			}
		
		}
		
		/**
		 * Returns results for all the column there were and are present in the collection
		 * 
		 * @return results
		 */
		
		List<Double> getSchneider() {
			
			assert this.nrOfUpdates == this.changesOfColumns && this.schneiderInitStat == true;
			
			return schneider;
			
		}
		
		/**
		 * Calculates Shenkin for all columns in the collection.
		 */
		
		void calculateShenkin() {
			
			
			for (int i = 0; i < cols.length; i++) {
				
				shenkin.add(cols[i].schneiderScore());
			}
		
		}
		
		/**
		 * Returns results for all the column there were and are present in the collection
		 * 
		 * @return results
		 */
		
		List<Double> getShenkin() {
			
			assert this.nrOfUpdates == this.changesOfColumns && this.shenkinInitStat == true;
			
			return shenkin;
			
		}

		/**
		 * Calculates Gerstein for all columns in the collection.
		 */
		
		void calculateGerstein() {
			
			
			for (int i = 0; i < cols.length; i++) {
				
				gerstein.add(cols[i].gersteinScore());
			}
		}
		
		/**
		 * Returns results for all the column there were and are present in the collection
		 * 
		 * @return results
		 */
		
		List<Double> getGerstein() {
			
			assert this.nrOfUpdates == this.changesOfColumns && this.gersteinInitStat == true;
			
			return gerstein;
			
		}
		
		/**
		 * Calculates TaylorNoGaps for all columns in the collection.
		 */

		void calculateTaylorNoGaps() {
			
			for (int i = 0; i < cols.length; i++) {
				
				taylorNoGaps.add(cols[i].SmallestTaylorSetNoGaps());
			}
		
		}
		
		
		/**
		 * Returns results for all the column there were and are present in the collection
		 * 
		 * @return results
		 */
		
		List<Integer> getTaylorNoGaps() {
			
			assert this.nrOfUpdates == this.changesOfColumns && this.taylorNoGapsInitSatat == true;
			
			return taylorNoGaps;
			
		}
		
		/**
		 * Calculates TayloGaps for all columns in the collection.
		 */
		
		void calculateTaylorGaps() {
			
			for (int i = 0; i < cols.length; i++) {
				
				taylorGaps.add(cols[i].SmallestTaylorSetGaps());
			}
		
		}
		
		/**
		 * Returns results for all the column there were and are present in the collection
		 * 
		 * @return results
		 */
		
		List<Integer> getTaylorGaps() {
			
			assert this.nrOfUpdates == this.changesOfColumns && this.taylorGapsInitStat == true;
			
			return taylorGaps;
			
		}

		/**
		 * Calculates Zvelibil for all columns in the collection.
		 */
		
		void calculateZvelibil() {
			
			for (int i = 0; i < cols.length; i++) {
				
				zvelibil.add(cols[i].zvelibilScore());
			}
		
		}
		
		/**
		 * Returns results for all the column there were and are present in the collection
		 * 
		 * @return results
		 */
		
		List<Integer> getZvelibil() {
			
			assert this.nrOfUpdates == this.changesOfColumns && this.zvelibilInitStat == true;
			
			return zvelibil;
			
		}
		
		/**
		 * Calculates Karlin for all columns in the collection.
		 */

		void calculateKarlin() {
			
			for (int i = 0; i < cols.length; i++) {
				
				karlin.add(cols[i].karlinScore());
			}
		
		}
		
		/**
		 * Returns results for all the column there were and are present in the collection
		 * 
		 * @return results
		 */
		
		List<Double> getKarlin() {
			
			assert this.nrOfUpdates == this.changesOfColumns && this.karlinInitStat == true;
			
			return karlin;
			
		}

		/**
		 * Calculates Armon for all columns in the collection.
		 */
		
		void calculateArmon() {
			
			for (int i = 0; i < cols.length; i++) {
				
				armon.add(cols[i].armonScore());
			}
		
		}
		
		/**
		 * Returns results for all the column there were and are present in the collection
		 * 
		 * @return results
		 */
		
		List<Double> getArmon() {
			
			assert this.nrOfUpdates == this.changesOfColumns && this.armonInitStat == true;
			
			return armon;
			
		}
		
		/**
		 * Calculates Thompson for all columns in the collection.
		 */

		void calculateThompson() {
			
			for (int i = 0; i < cols.length; i++) {
				
				thompson.add(cols[i].thompsonScore());
			}
		
		}
		
		/**
		 * Returns results for all the column there were and are present in the collection
		 * 
		 * @return results
		 */
		
		List<Double> getThompson() {
			
			assert this.nrOfUpdates == this.changesOfColumns && this.thompsonInitStat == true;
			
			return thompson;
			
		}
		
		/**
		 * Calculates Lancet for all columns in the collection.
		 */

		void calculateLancet() {
			
			for (int i = 0; i < cols.length; i++) {
				
				lancet.add(cols[i].lancetScore());
			}
		
		}
		
		/**
		 * Returns results for all the column there were and are present in the collection
		 * 
		 * @return results
		 */
		
		List<Double> getLancet() {
	
			assert this.nrOfUpdates == this.changesOfColumns && this.lancetInitStat == true;
			
			return lancet;
			
		}
		
		/**
		 * Calculates Mirny for all columns in the collection.
		 */

		void calculateMirny() {
			
			for (int i = 0; i < cols.length; i++) {
				
				mirny.add(cols[i].mirnyScore());
			}
		
		}
		
		/**
		 * Returns results for all the column there were and are present in the collection
		 * 
		 * @return results
		 */
		
		List<Double> getMirny() {
			
			assert this.nrOfUpdates == this.changesOfColumns && this.mirnyInitStat == true; 
			
			return mirny;
			
		}
		
		/**
		 * Calculates Wiliamson for all columns in the collection.
		 */

		void calculateWilliamson() {
			
			for (int i = 0; i < cols.length; i++) {
				
				williamson.add(cols[i].williamsonScore());
			}
		
		}
		
		/**
		 * Returns results for all the column there were and are present in the collection
		 * 
		 * @return results
		 */
		
		List<Double> getWilliamson() {
			
			assert this.nrOfUpdates == this.changesOfColumns && this.williamsonInitStat == true;
			
			return williamson;
			
		}

		/**
		 * Calculates Landgraf for all columns in the collection.
		 */
		
		void calculateLandgarf() {
			
			for (int i = 0; i < cols.length; i++) {
				
				landgraf.add( cols[i].landgrafScore());
			}
		
		}
		
		/**
		 * Returns results for all the column there were and are present in the collection
		 * 
		 * @return results
		 */
		
		List<Double> getLandgraf() {
			
			assert this.nrOfUpdates == this.changesOfColumns && this.landgrafInitStat == true;
			
			return landgraf;
			
		}
		
		/**
		 * Calculates Sander for all columns in the collection.
		 */

		void calculateSander() {
			
			for (int i = 0; i < cols.length; i++) {
				
				sander.add(cols[i].sanderScore());
			}
		
		}
		
		/**
		 * Returns results for all the column there were and are present in the collection
		 * 
		 * @return results
		 */
		
		List<Double> getSander() {
			
			assert this.nrOfUpdates == this.changesOfColumns && this.sanderInitStat == true;
			
			return sander;
			
		}
		
		/**
		 * Calculates Valdar for all columns in the collection.
		 */
		
		void calculateValdar() {
			
			for (int i = 0; i < cols.length; i++) {
				
				valdar.add(cols[i].valdarScore());
			}
		
		}
		
		/**
		 * Returns results for all the column there were and are present in the collection
		 * 
		 * @return results
		 */
		
		List<Double> getValdar() {
			
			assert this.nrOfUpdates == this.changesOfColumns && this.valdarInitStat == true;
			
			return valdar;
			
		}

		/**
		 * Calculates Jores for all columns in the collection.
		 */
		
		void calculateJores() {
			
			for (int i = 0; i < cols.length; i++) {
				
				jores.add(cols[i].joresScore());
				
			}
			
		}
		
		/**
		 * Returns results for all the column there were and are present in the collection
		 * 
		 * @return results
		 */
		
		List<Double> getJores() {
			
			assert this.nrOfUpdates == this.changesOfColumns && this.joresInitStat == true;
			
			return jores;
		}
		
		/**
		 * Changes columns in the collection. Replaces current set of columns with a new one.
		 * @param changeStart
		 */
		
		void changeColumns(int changeStart) {
			
			assert changeStart < this.matrix.numberOfColumns();
			
			assert matrix.numberOfColumns() - ((this.changesOfColumns + 1) * this.cols.length) > this.cols.length;
			
			for (int i = 0; i < cols.length; i++) {
				
				cols[i] = new Column(matrix, changeStart + i);
			}
			
			changesOfColumns++;
			
		}
		
		/**
		 * Updates the calculations in the collection after columns have been changed.
		 * Has to be performed after each change.
		 * LOOK OUT: Needs exactly the same parameter values as calculationInitialization().
		 *  
		 * @param kabat
		 * @param jores
		 * @param schneider
		 * @param shenkin
		 * @param gerstein
		 * @param taylorNoGaps
		 * @param taylorGaps
		 * @param zvelibil
		 * @param karlin
		 * @param armon
		 * @param thompson
		 * @param lancet
		 * @param mirny
		 * @param williamson
		 * @param landgraf
		 * @param sander
		 * @param valdar
		 */
		
		void calculationUpdate(boolean kabat, boolean jores, boolean schneider, boolean shenkin, boolean gerstein, boolean taylorNoGaps, boolean taylorGaps, boolean zvelibil, boolean karlin, boolean armon, boolean thompson, boolean lancet, boolean mirny, boolean williamson, boolean landgraf, boolean sander, boolean valdar) {
			
			assert this.calcsInitialized == true;
			
			assert this.kabatInitStat == kabat && this.joresInitStat == jores && this.schneiderInitStat == schneider && this.shenkinInitStat == shenkin && this.gersteinInitStat == gerstein && this.taylorNoGapsInitSatat == taylorNoGaps && this.taylorGapsInitStat == taylorGaps && this.zvelibilInitStat == zvelibil && this.karlinInitStat == karlin && this.armonInitStat == armon && this.thompsonInitStat == thompson && this.lancetInitStat == lancet && this.mirnyInitStat == mirny && this.williamsonInitStat == williamson && this.landgrafInitStat == landgraf && this.sanderInitStat == sander && this.valdarInitStat == valdar;
			
			this.nrOfUpdates++;
			
			if (kabat) {
				
				this.calculateKabat();
				
				
				}
			
			if(jores) {
				
				this.calculateJores();
				
			}
			
			if(schneider) {
				
				this.calculateSchneider();
				
			}
			
			if(shenkin) {
				
				this.calculateShenkin();
				
			}
			
			if(gerstein) {
				
				this.calculateGerstein();
				
			}
			
			if (taylorNoGaps) {
				
				this.calculateTaylorNoGaps();
			
			}
			
			if (taylorGaps) {
				
				this.calculateTaylorGaps();
				
			}
			
			if (zvelibil) {
				
				this.calculateZvelibil();
				
			}
			
			if (karlin) {
				
				this.calculateKarlin();
				
			}
			
			if(armon) {
				
				this.calculateArmon();
				
			}
			
			if(thompson) {
				
				this.calculateThompson();
				
			}
			
			if(lancet) {
				
				this.calculateLancet();
				
			}
			
			if(mirny) {
				
				this.calculateMirny();
				
			}
			
			if(williamson) {
				
				this.calculateWilliamson();
				
			}
			
			if(landgraf) {
				
				this.calculateLandgarf();
				
			}
			
			if(sander) {
				
				this.calculateSander();
				
				
			}
			
			if(valdar) {
				
				this.calculateValdar();
			
				
			}
			
		}
		
		/**
		 * Has to be performed after object formation if any calculations are to be performed.
		 * Can be performed only once per object creation.
		 * 
		 * @param kabat
		 * @param jores
		 * @param schneider
		 * @param shenkin
		 * @param gerstein
		 * @param taylorNoGaps
		 * @param taylorGaps
		 * @param zvelibil
		 * @param karlin
		 * @param armon
		 * @param thompson
		 * @param lancet
		 * @param mirny
		 * @param williamson
		 * @param landgraf
		 * @param sander
		 * @param valdar
		 */
		
		void calculationInitializtion(boolean kabat, boolean jores, boolean schneider, boolean shenkin, boolean gerstein, boolean taylorNoGaps, boolean taylorGaps, boolean zvelibil, boolean karlin, boolean armon, boolean thompson, boolean lancet, boolean mirny, boolean williamson, boolean landgraf, boolean sander, boolean valdar) {
			
			assert !(this.calcsInitialized == true);
			
			if (kabat) {
				
				this.kabat = new ArrayList<Double>();
				
				this.calculateKabat();
				
				this.kabatInitStat = kabat;
				
				}
			
			if(jores) {
				
				this.jores = new ArrayList<Double>();
				
				this.calculateJores();
				
				this.joresInitStat = jores;
				
			}
			
			if(schneider) {
				
				this.schneider = new ArrayList<Double>();
				
				this.calculateSchneider();
				
				this.schneiderInitStat = schneider;
			
			}
			
			if(shenkin) {
				
				this.shenkin = new ArrayList<Double>();
				
				this.calculateShenkin();
				
				this.shenkinInitStat = shenkin;
				
			}
			
			if(gerstein) {
				
				this.gerstein = new ArrayList<Double>();
				
				this.calculateGerstein();
				
				this.gersteinInitStat = gerstein;
				
			}
			
			if (taylorNoGaps) {
				
				this.taylorNoGaps = new ArrayList<Integer>();
				
				this.calculateTaylorNoGaps();
				
				this.taylorNoGapsInitSatat = true;
			}
			
			if (taylorGaps) {
				
				this.taylorGaps = new ArrayList<Integer>();
				
				this.calculateTaylorGaps();
				
				this.taylorGapsInitStat = taylorGaps;
				
			}
			
			if (zvelibil) {
				
				this.zvelibil = new ArrayList<Integer>();
				
				this.calculateZvelibil();
				
				this.zvelibilInitStat = zvelibil;
				
			}
			
			if (karlin) {
				
				this.karlin = new ArrayList<Double>();
				
				this.calculateKarlin();
				
				this.karlinInitStat = karlin;
				
			}
			
			if(armon) {
				
				this.armon = new ArrayList<Double>();
				
				this.calculateArmon();
				
				this.armonInitStat = armon;
			
			}
			
			if(thompson) {
				
				this.thompson = new ArrayList<Double>();
				
				this.calculateThompson();
				
				this.thompsonInitStat = thompson;
				
			}
			
			if(lancet) {
				
				this.lancet = new ArrayList<Double>();
				
				this.calculateLancet();
				
				this.lancetInitStat = lancet;
				
			}
			
			if(mirny) {
				
				this.mirny = new ArrayList<Double>();
				
				this.calculateMirny();
				
				this.mirnyInitStat = mirny;
			
			}
			
			if(williamson) {
				
				this.williamson = new ArrayList<Double>();
				
				this.calculateWilliamson();
				
				this.williamsonInitStat = williamson;
				
			}
			
			if(landgraf) {
				
				this.landgraf = new ArrayList<Double>();
				
				this.calculateLandgarf();
				
				this.landgrafInitStat = landgraf;
				
			}
			
			if(sander) {
				
				this.sander = new ArrayList<Double>();
				
				this.calculateSander();
				
				this.sanderInitStat = sander;
				
				
			}
			
			if(valdar) {
				
				this.valdar = new ArrayList<Double>();
				
				this.calculateValdar();
				
				this.valdarInitStat = valdar;
		
				
			}
		
			this.calcsInitialized = true;
		}
		
		
		
		
	}




