package compbio.conservation;

import java.util.ArrayList;
import java.util.List;
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

public class _ColumnCollectionShortList {
	
			
			private int changesOfColumns = 0;
		
			/**
			 * Keeps track of the number of   updates
			 */
			
			private int nrOfUpdates = 0;
			/**
			 * Holds the columns currently present in the window.
			 */
			private final List<_Column> cols;
			
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

			
			public _ColumnCollectionShortList(AminoAcidMatrix m, int nrOfColumns) {
			
			if (m == null) {
				
				throw new IllegalArgumentException("m must not be null");
			}
			
			assert m.numberOfColumns() != 0 : "Sth wrong with the matrix, not null but doesn't have any residues";
			
			assert nrOfColumns < m.numberOfColumns();
			
			matrix = m;
			
			cols = new ArrayList<_Column>(nrOfColumns);
			
			for (int i = 0; i < nrOfColumns; i++) {
				
				cols.add(new _Column(matrix, i));
				
			}
			
			}	
			
			/**
			 * Returns columns currently present in the collection.
			 * 
			 * @return column collection
			 */
			
			List<_Column> getColumnCollection () {
			
			return cols;
			
			}
			
			/**
			 * Returns collection lenght( window width).
			 * 
			 * @return collection length
			 */

			int collectionLength() {
			
			int len = cols.size();
			
			return len;
			
			}
			
			/**
			 * Calculates Kabat for all columns in the collection.
			 */

			void calculateKabat(int nrOfUpdates) {
				
				if (nrOfUpdates == 0) {
				
				for (int i = 0; i < cols.size(); i++) {
					
					kabat.add(cols.get(i).kabatScore());
				}
				
				}
				
				else {
					
					kabat.remove(0);
					
					kabat.add(cols.get(this.cols.size() - 1).kabatScore());
				}
			
			}
			
			/**
			 * Returns results for all the current columns in the collection
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
			
			void calculateSchneider(int nrOfUpdates) {
				
				if (nrOfUpdates == 0) {
				
				for (int i = 0; i < cols.size(); i++) {
					
					schneider.add(cols.get(i).schneiderScore());
				}
				
				}
				
				else {
					
					schneider.remove(0);
					
					schneider.add(cols.get(this.cols.size() - 1).schneiderScore());
					
				}
			
			}
			
			/**
			 * Returns results for all the current columns in the collection
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
			
			void calculateShenkin(int nrOfUpdates) {
				
				if (nrOfUpdates == 0) {
				
				for (int i = 0; i < cols.size(); i++) {
					
					shenkin.add(cols.get(i).shenkinScore());
				}
				
				}
				
				else {
					
					this.shenkin.remove(0);
					
					this.shenkin.add(cols.get(this.cols.size() - 1).shenkinScore());
				}
			
			}
			
			/**
			 * Returns results for all the current columns in the collection
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
			
			void calculateGerstein(int nrOfUpdaes ) {
				
				if (nrOfUpdates == 0) {
				
				for (int i = 0; i < cols.size(); i++) {
					
					gerstein.add(cols.get(i).gersteinScore());
				}
				
				}
				
				else {
					
					gerstein.remove(0);
					
					gerstein.add(cols.get(this.cols.size() - 1).gersteinScore());
				}
			}
			
			/**
			 * Returns results for all the current columns in the collection
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

			void calculateTaylorNoGaps(int nrOfUpdates) {
				
				if (nrOfUpdates == 0) {
				
				for (int i = 0; i < cols.size(); i++) {
					
					taylorNoGaps.add(cols.get(i).SmallestTaylorSetNoGaps());
				}
				
				}
				
				else {
					
					this.taylorNoGaps.remove(0);
					
					this.taylorNoGaps.add(cols.get(this.cols.size() - 1).SmallestTaylorSetNoGaps());
				}
			
			}
			
			
			/**
			 * Returns results for all the current columns in the collection
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
			
			void calculateTaylorGaps(int nrOfUpdates) {
				
				if(nrOfUpdates == 0) {
				
				for (int i = 0; i < cols.size(); i++) {
					
					taylorGaps.add(cols.get(i).SmallestTaylorSetGaps());
				}
				
				}
				
				else {
					
					this.taylorNoGaps.remove(0);
					
					this.taylorNoGaps.add(cols.get(this.cols.size() - 1).SmallestTaylorSetNoGaps());
				}
			
			}
			
			/**
			 * Returns results for all the current columns in the collection
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
			
			void calculateZvelibil(int nrOfUpdates) {
				
				if (nrOfUpdates == 0){
				
				for (int i = 0; i < cols.size(); i++) {
					
					zvelibil.add(cols.get(i).zvelibilScore());
				}
				
				}
				
				else {
					
					this.zvelibil.remove(0);
					
					this.zvelibil.add(cols.get(this.cols.size() - 1).zvelibilScore());
				}
			
			}
			
			/**
			 * Returns results for all the current columns in the collection
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

			void calculateKarlin(int nrOfUpdates) {
				
				if(nrOfUpdates == 0) {
				
				for (int i = 0; i < cols.size(); i++) {
					
					karlin.add(cols.get(i).karlinScore());
				}
				
				}
				
				else {
					
					this.karlin.remove(0);
					
					this.karlin.add(this.cols.get(this.cols.size() - 1).karlinScore());
					
				}
				}
			
			/**
			 * Returns results for all the current columns in the collection
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
			
			void calculateArmon(int nrOfUpdates) {
				
				if (nrOfUpdates == 0){
				
				for (int i = 0; i < cols.size(); i++) {
					
					armon.add(cols.get(i).armonScore());
				}
				
				}
				
				else {
					
					this.armon.remove(0);
					
					this.armon.add(this.cols.get(this.cols.size() - 1).armonScore());
					
				}
			
			}
			
			/**
			 * Returns results for all the current columns in the collection
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

			void calculateThompson(int nrOfUdates) {
				
				if(nrOfUpdates == 0) {
				
				for (int i = 0; i < cols.size(); i++) {
					
					thompson.add(cols.get(i).thompsonScore());
				}
				
				}
				
				else {
					
					this.thompson.remove(0);
					
					this.thompson.add(this.cols.get(this.cols.size() - 1).thompsonScore());
				}
			
			}
			
			/**
			 * Returns results for all the current columns in the collection
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

			void calculateLancet(int nrOfUpdates) {
				
				if (nrOfUpdates == 0) {
				
				for (int i = 0; i < cols.size(); i++) {
					
					lancet.add(cols.get(i).notLancetScore());
				}
				
				}
			
				else {
					
					this.lancet.remove(0);
					
					this.lancet.add(this.cols.get(this.cols.size() - 1).notLancetScore());
				}
			}
			
			/**
			 * Returns results for all the current columns in the collection
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

			void calculateMirny(int nrOfUpdates) {
				
				if (nrOfUpdates == 0) {
				
				for (int i = 0; i < cols.size(); i++) {
					
					mirny.add(cols.get(i).mirnyScore());
				}
				
				}
				
				else {
					
					this.mirny.remove(0);
					
					this.mirny.add(this.cols.get(this.cols.size() - 1).mirnyScore());
					
				}
			
			}
			
			/**
			 * Returns results for all the current columns in the collection
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

			void calculateWilliamson(int nrOfUdates) {
				
				if(nrOfUpdates == 0) {
				
				for (int i = 0; i < cols.size(); i++) {
					
					williamson.add(cols.get(i).williamsonScore());
				}
				
				}
				
				else {
					
					this.williamson.remove(0);
					
					this.williamson.add(this.cols.get(this.cols.size() - 1).williamsonScore());
				}
			
			}
			
			/**
			 * Returns results for all the current columns in the collection
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
			
			void calculateLandgarf(int nrOfUpdates) {
				
				if (nrOfUpdates == 0){
				
				for (int i = 0; i < cols.size(); i++) {
					
					landgraf.add( cols.get(i).landgrafScore());
				}
				
				}
				
				else {
					
					this.landgraf.remove(0);
					
					this.landgraf.add(this.cols.get(this.cols.size() - 1).landgrafScore());
				}
			
			}
			
			/**
			 * Returns results for all the current columns in the collection
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

			void calculateSander(int nrOfUpdates) {
				
				if(nrOfUpdates == 0) {
				
				for (int i = 0; i < cols.size(); i++) {
					
					sander.add(cols.get(i).sanderScore());
				}
				
				}
				
				else {
					
					this.sander.remove(0);
					
					this.sander.add(this.cols.get(this.cols.size() - 1).sanderScore());
				}
			
			}
			
			/**
			 * Returns results for all the current columns in the collection
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
			
			void calculateValdar(int nrOfUpdates) {
				
				if (nrOfUpdates == 0) {
				
				for (int i = 0; i < cols.size(); i++) {
					
					valdar.add(cols.get(i).valdarScore());
				}
				
				}
				
				else {
					
					this.valdar.remove(0);
					
					this.valdar.add(this.cols.get(this.cols.size() -1).valdarScore());
				}
			
			}
			
			/**
			 * Returns results for all the current columns in the collection
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
			
			void calculateJores(int nrOfUpdates) {
				
				if (nrOfUpdates == 0) {
				
				for (int i = 0; i < cols.size(); i++) {
					
					jores.add(cols.get(i).joresScore());
					
				}
				
				}
				
				else {
					
					this.jores.remove(0);
					
					this.jores.add(this.cols.get(this.cols.size() - 1).joresScore());
				}
				
			}
			
			/**
			 * Returns results for all the current columns in the collection
			 * 
			 * @return results
			 */
			
			List<Double> getJores() {
				
				assert this.nrOfUpdates == this.changesOfColumns && this.joresInitStat == true;
				
				return jores;
			}
			
			/**
			 * Changes columns in the collection. Replaces current set of columns with a new one.
			 * 
			 * @param newColumnIndex index of the new column which replaces the oldest in the collection
			 */
			
			void changeColumns(int newColumnIndex) {
				
				assert newColumnIndex < this.matrix.numberOfColumns();
				
				//assert matrix.numberOfColumns() - ((this.changesOfColumns + 1) * this.cols.length) > this.cols.length;
				
					this.cols.remove(0);
				
					cols.add(new _Column(matrix, newColumnIndex));
				
				changesOfColumns++;
				
				
			}
			
			/**
			 * Updates the calculations in the collection after columns have been changed.
			 * Therefore called only after changeOfColumns() method.
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
					
					this.calculateKabat(this.changesOfColumns);
					
					
					}
				
				if(jores) {
					
					this.calculateJores(this.changesOfColumns);
					
				}
				
				if(schneider) {
					
					this.calculateSchneider(this.changesOfColumns);
					
				}
				
				if(shenkin) {
					
					this.calculateShenkin(this.changesOfColumns);
					
				}
				
				if(gerstein) {
					
					this.calculateGerstein(this.changesOfColumns);
					
				}
				
				if (taylorNoGaps) {
					
					this.calculateTaylorNoGaps(this.changesOfColumns);
				
				}
				
				if (taylorGaps) {
					
					this.calculateTaylorGaps(this.changesOfColumns);
					
				}
				
				if (zvelibil) {
					
					this.calculateZvelibil(this.changesOfColumns);
					
				}
				
				if (karlin) {
					
					this.calculateKarlin(this.changesOfColumns);
					
				}
				
				if(armon) {
					
					this.calculateArmon(this.changesOfColumns);
					
				}
				
				if(thompson) {
					
					this.calculateThompson(this.changesOfColumns);
					
				}
				
				if(lancet) {
					
					this.calculateLancet(this.changesOfColumns);
					
				}
				
				if(mirny) {
					
					this.calculateMirny(this.changesOfColumns);
					
				}
				
				if(williamson) {
					
					this.calculateWilliamson(this.changesOfColumns);
					
				}
				
				if(landgraf) {
					
					this.calculateLandgarf(this.changesOfColumns);
					
				}
				
				if(sander) {
					
					this.calculateSander(this.changesOfColumns);
					
					
				}
				
				if(valdar) {
					
					this.calculateValdar(this.changesOfColumns);
				
					
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
					
					this.kabat = new ArrayList<Double>(this.cols.size());
					
					this.calculateKabat(this.changesOfColumns);
					
					this.kabatInitStat = kabat;
					
					}
				
				if(jores) {
					
					this.jores = new ArrayList<Double>(this.cols.size());
					
					this.calculateJores(this.changesOfColumns);
					
					this.joresInitStat = jores;
					
				}
				
				if(schneider) {
					
					this.schneider = new ArrayList<Double>(this.cols.size());
					
					this.calculateSchneider(this.changesOfColumns);
					
					this.schneiderInitStat = schneider;
				
				}
				
				if(shenkin) {
					
					this.shenkin = new ArrayList<Double>(this.cols.size());
					
					this.calculateShenkin(this.changesOfColumns);
					
					this.shenkinInitStat = shenkin;
					
				}
				
				if(gerstein) {
					
					this.gerstein = new ArrayList<Double>(this.cols.size());
					
					this.calculateGerstein(this.changesOfColumns);
					
					this.gersteinInitStat = gerstein;
					
				}
				
				if (taylorNoGaps) {
					
					this.taylorNoGaps = new ArrayList<Integer>(this.cols.size());
					
					this.calculateTaylorNoGaps(this.changesOfColumns);
					
					this.taylorNoGapsInitSatat = taylorNoGaps;
				}
				
				if (taylorGaps) {
					
					this.taylorGaps = new ArrayList<Integer>(this.cols.size());
					
					this.calculateTaylorGaps(this.changesOfColumns);
					
					this.taylorGapsInitStat = taylorGaps;
					
				}
				
				if (zvelibil) {
					
					this.zvelibil = new ArrayList<Integer>(this.cols.size());
					
					this.calculateZvelibil(this.changesOfColumns);
					
					this.zvelibilInitStat = zvelibil;
					
				}
				
				if (karlin) {
					
					this.karlin = new ArrayList<Double>(this.cols.size());
					
					this.calculateKarlin(this.changesOfColumns);
					
					this.karlinInitStat = karlin;
					
				}
				
				if(armon) {
					
					this.armon = new ArrayList<Double>(this.cols.size());
					
					this.calculateArmon(this.changesOfColumns);
					
					this.armonInitStat = armon;
				
				}
				
				if(thompson) {
					
					this.thompson = new ArrayList<Double>(this.cols.size());
					
					this.calculateThompson(this.changesOfColumns);
					
					this.thompsonInitStat = thompson;
					
				}
				
				if(lancet) {
					
					this.lancet = new ArrayList<Double>(this.cols.size());
					
					this.calculateLancet(this.changesOfColumns);
					
					this.lancetInitStat = lancet;
					
				}
				
				if(mirny) {
					
					this.mirny = new ArrayList<Double>(this.cols.size());
					
					this.calculateMirny(this.changesOfColumns);
					
					this.mirnyInitStat = mirny;
				
				}
				
				if(williamson) {
					
					this.williamson = new ArrayList<Double>(this.cols.size());
					
					this.calculateWilliamson(this.changesOfColumns);
					
					this.williamsonInitStat = williamson;
					
				}
				
				if(landgraf) {
					
					this.landgraf = new ArrayList<Double>(this.cols.size());
					
					this.calculateLandgarf(this.changesOfColumns);
					
					this.landgrafInitStat = landgraf;
					
				}
				
				if(sander) {
					
					this.sander = new ArrayList<Double>(this.cols.size());
					
					this.calculateSander(this.changesOfColumns);
					
					this.sanderInitStat = sander;
					
					
				}
				
				if(valdar) {
					
					this.valdar = new ArrayList<Double>(this.cols.size());
					
					this.calculateValdar(this.changesOfColumns);
					
					this.valdarInitStat = valdar;
			
					
				}
			
				this.calcsInitialized = true;
			}
			
			_Column getNewestColumn() {
				
				return this.cols.get(this.cols.size() - 1);
			}
		
			
}

	





