package compbio.conservation;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class _ColumnWithoutCollection {

	/**
	 * Array consisting of character representation of amino acids in the column
	 */
	private final char[] columnArr;
	
	/**
	 * Maps a characters representing amino acids with their occurrence.
	 */
    private final Map<Character,Integer> acidsIntMap;

    /** 
     * Reference to the matrix object the column was taken from. Indices of the amino acids in the column correspond to the indices of the the sequences they belong to in this matrix.
     */
    private final AminoAcidMatrix matrix;
    
    /**
	 * Holds the of Kabat result.
	 */
	
	private double kabat = 0;
	
	/**
	 * Kabat initialization status
	 */
	
	private boolean kabatInitStat = false;
	
	/**
	 * Holds Jores result.
	 */
	
	private double jores = 0;
	
	/**
	 * Jores initialization status.
	 */
	
	private boolean joresInitStat = false;
	
	/**
	 * Holds Schneider result.
	 */

	private double schneider = 0;
	
	/**
	 * Schneider initialization status.
	 */
	
	private boolean schneiderInitStat = false;
	
	/**
	 * Holds Shenkin result.
	 */

	private double shenkin = 0;
	
	/**
	 * Shenkin initialization status.
	 */
	
	private boolean shenkinInitStat = false;
	
	/**
	 * Holds Gerstein result.
	 */

	private double gerstein = 0;
	
	/**
	 * Gerstein initialization status.
	 */
	
	private boolean gersteinInitStat = false;
	
	/**
	 * Holds Taylor result.
	 */

	private int taylorNoGaps = 0;
	
	/**
	 * Taylor initialization status.
	 */
	
	private boolean taylorNoGapsInitSatat = false;
	
	/**
	 * Holds Taylor result.
	 */

	private int taylorGaps = 0;
	
	/**
	 * Taylor initialization status.
	 */
	
	private boolean taylorGapsInitStat = false;
	
	/**
	 * Holds Zvelibil result.
	 */

	private int zvelibil = 0;
	
	/**
	 * Zvelibil initialization status.
	 */
	
	private boolean zvelibilInitStat = false;
	
	/**
	 * Holds Karlin result.
	 */

	private double karlin = 0;
	
	/** 
	 * Karlin initialization status.
	 */
	
	private boolean karlinInitStat = false;
	
	/**
	 * Holds Armon result.
	 */

	private double armon = 0;
	
	/**
	 * Armon initialization status.
	 */
	
	private boolean armonInitStat = false;
	
	/**
	 * Holds Thompson Result.
	 */

	private double thompson = 0;
	
	/**
	 * Thompson initialization status.
	 */
	
	private boolean thompsonInitStat = false;
	
	/**
	 * Holds Lancet result.
	 */

	private double lancet = 0;
	
	/**
	 * Lancet initialization status.
	 */
	
	private boolean lancetInitStat = false;
	
	/**
	 * Holds Mirny result.
	 */

	private double mirny = 0;
	
	/**
	 * Mirny initialization status.
	 */
	
	private boolean mirnyInitStat = false;
	
	/**
	 * List Williamson result.
	 */

	private double williamson = 0;
	
	/**
	 * Williamson initialization status.
	 */
	
	private boolean williamsonInitStat = false;
	
	/**
	 * Holds landgraf result.
	 */

	private double landgraf = 0;
	
	/**
	 * Landgraf initialization status.
	 */
	
	private boolean landgrafInitStat = false;
	
	/**
	 * Holds Sander result.
	 */

	private double sander = 0;
	
	/**
	 * Sander initialization status.
	 */
	
	private boolean sanderInitStat = false;
	
	/**
	 * Holds  Valdar result.
	 */

	private double valdar = 0;
	
	/**
	 * Valdar initialization status.
	 */
	
	private boolean valdarInitStat = false;
	
    /**
     * Constructor that allows manual creation of the column of the length of 10.
     * @param a
     * @param b
     * @param c
     * @param d
     * @param e
     * @param f
     * @param g
     * @param h
     * @param i
     * @param j
     * @param m
     */

    public _ColumnWithoutCollection(char a, char b, char c, char d, char e, char f,
	  char g, char h, char i, char j, AminoAcidMatrix m) {

	columnArr = new char[10];

	columnArr[0] = a;
	columnArr[1] = b;
	columnArr[2] = c;
	columnArr[3] = d;
	columnArr[4] = e;
	columnArr[5] = f;
	columnArr[6] = g;
	columnArr[7] = h;
	columnArr[8] = i;
	columnArr[9] = j;


    acidsIntMap = Alphabet.calculateOccurance(columnArr);
    
    assert !acidsIntMap.isEmpty() && acidsIntMap != null;
    
    matrix = m;
    
    }


    /**
     * Constructor. Creates a column object by copying the appropriate column form the matrix.
     * Takes a reference to the alignment matrix and the index of the corresponding column in the matrix. 
     * @param mmatrix
     * @param columnNr
     */
    public _ColumnWithoutCollection(AminoAcidMatrix matrix, int columnNr) {
    	 

	if (matrix == null) {
	    throw new IllegalArgumentException("Matrix must not be null");
	}
	
	assert columnNr >= 0 && columnNr <  matrix.numberOfColumns();
	
	this.matrix = matrix;

	String columnStr = new String(matrix.getColumn(columnNr));

	String upColumnString = columnStr.toUpperCase();

	columnArr = upColumnString.toCharArray();
	
	Map<Character, Integer> map = Alphabet.calculateOccurance(columnArr);
	
	acidsIntMap = Collections.unmodifiableMap(map);
	
	assert acidsIntMap != null && !acidsIntMap.isEmpty();
	
	}
	/**
	 * Checks if all but one residues in the column are gaps.
	 * 
	 * @return true if all but one residues are gaps, false if not 
	 */
	
    boolean allButOneGaps() {
    	
    	if (acidsIntMap.containsKey('-') && acidsIntMap.get('-') == columnArr.length - 1) {
			
			return true;
		}
	
	else {
		
		return false;
		
	}
	
    }
    /**
     * Checks whether there is only one residue type in the column. Gap is not counted as a residue type.
     * 
     * @return returns true if there is only one residue type in the column, false if not
     */
    boolean oneResidueTypeNoGaps() {
    	
    	if (acidsIntMap.size() == 1 && acidsIntMap.containsKey('-') == false) {
    		
    		return true;
    		
    	}
    	
    	else {
    		
    		return false;
    		
    	}
    	
    }
    
    
    /** 
     * Checks whether column contains gaps.
     * 
     * @return true if column contains gaps
     */
    
    boolean containsGaps() {
    	
    	if (acidsIntMap.containsKey('-')) {
    		
    		return true;
    	}
    	
    	else {
    		
    		return false;
    		
    	}
    	
    }
 
    /**
     * Counts the number of different amino acids in the column.
     * Gap is counted as a 21 amino acid.
     * 
     * @return the number of different amino acids
     */
    
    int numberOfAcidsWithGap() {
    	
    	return acidsIntMap.size(); 
    	
    }
    
    /**
     * Counts the number of different amino acids in the column.
     * Gap not counted as a 21 amino acid.
     * 
     * @return the number of different amino acids
     */
    
    int numberOfAcidsNoGap() {
    
    if(this.containsGaps() == true) {
    	
    	return acidsIntMap.size() - 1;
    }
    
    else {
    	
    	return acidsIntMap.size();
    
    }
    
    }
    /** 
     * Calculates the number of times the most common amino acid in the column occurs.
     * Does not count gap as 21 amino acid.
     * 
     * @return times the most common amino acid occurs
     */
    
    int mostCommonNumber() {
    	
    	int max = 0;
    	
    	Set<Character> keys = acidsIntMap.keySet();
    	
    	Iterator<Character> itr = keys.iterator();
    	
    	while(itr.hasNext()) {
    		
    	Character key = itr.next();
    		
    		if (key != '-' && acidsIntMap.get(key) > max) {
    			
    			max = acidsIntMap.get(key);
    			
    		}
    	
    	}
    	
    	assert max != 0 : "Zero in the most Common Number";
    	
    	return max;
    	
    }
    
    Map<Character,Integer> getAcidsIntMap() {
    
    return acidsIntMap;
    
    
    }
    /**
     * Gets length of the column.
     * 
     * @return length of the column
     */
    			
    int length() {
	 
	 return columnArr.length;
	 
    }
 
    /**
     * Calculates Kabat score for the column.
     *  
     * @return Kabat score
     */
    void kabatScore() {
	 
        assert this.mostCommonNumber() > 0 && this.mostCommonNumber() < columnArr.length + 1;	 
         		
		double result = this.length() * this.numberOfAcidsNoGap()/ this.mostCommonNumber(); 
	
		this.kabat = result;
		
		this.kabatInitStat = true;
		
	}
    
    /**
     * Gets Kabat score.
     * 
     * @return Kabat score 
     */
    
    double getKabat() {
    	
    	if (this.kabatInitStat == false) {
    		
    		this.kabatScore();
    		
    		return this.kabat;
    	}
    	
    	else {
    		
    		return this.kabat;
    	}
    }
 
    /**
     * Calculates Jores score for the column.
     * When calculating the number of distinct pair formed calculates a combination of two out of the aminoacid types
     * present in the column. Than adds the self pairs. However when there is only one of a type
     * in the column ne self pair is formed. 
     * 
     * @return Jores score 
     */
    void joresScore() {
	 
	 double result = 0.0;
	 
	 // special case #1 one residue type only
	 
	 if (this.oneResidueTypeNoGaps() == true) {
		 
		 result = 1.0;
		 
		 this.jores = result;
	 
	 }
	
	 // special case #2 all but one are gaps
	 
	 if (this.allButOneGaps() == true) {
		 
		 this.jores = result;
		 
	 }
	 
	 int samePairs = 0;
	 
	 int differentPairs = 0;
	 
	 Map<Character,Integer> acidsIntMapCopy = new HashMap<Character,Integer>(acidsIntMap);
	 
	 Set<Character> keys = acidsIntMapCopy.keySet();
	 
	 keys.remove('-');
	 
	 int types = keys.size();
	 
	 differentPairs = types * (types - 1) / 2;
	 
	 Iterator<Integer> itr = acidsIntMapCopy.values().iterator();
	 
	 while (itr.hasNext()) {
		 
		 if (itr.next() > 1) {
			 
			 samePairs++;
			 
		 }
		 
	 }
	 
	 int totalPairs = samePairs + differentPairs;
	 
	 Iterator<Character> itr2 = acidsIntMapCopy.keySet().iterator();
	 
	 int max1 = 0;
	 
	 Character maxKey = null;
	 
	 while(itr2.hasNext()) {
		 
		 Character key = itr2.next();
		 
		 if (acidsIntMapCopy.get(key) > max1) {
			 
			 maxKey = key;
			 
			 max1 = acidsIntMapCopy.get(key);
			 
		 }
		 
	 }
	 
	 acidsIntMapCopy.remove(maxKey);
	 
	 Iterator<Integer> itr3 = acidsIntMapCopy.values().iterator();
	 
	 int max2 = 0;
	 
	 while(itr3.hasNext()) {
		 
		 int value = itr3.next();
		 
		 if(value > max2) {
			 
			 max2 = value;
			 
		 }
		 
	 }
	 
	 int mostFreqNr = 0;
	 
	 if (max2 == 0) {
		 
		 mostFreqNr = (max1) * (max1 - 1) / 2;
		 
	 }
	 
	 else {
		 
		 if (max1 == max2) {
			 
			 mostFreqNr = max1 * max2;
			 
		 }
		 
		 else {
			 
			 int same = (max1) * (max1 - 1) / 2;
			 
			 int diff = max1 * max2;
			 
			 if (same > diff) {
				 
				 mostFreqNr = same;
				 
			 }
				 
			 else {
					 
					 mostFreqNr = diff;
				 }
			 
		 	}
		 
		 this.joresInitStat = true;
		 
	 	}
		 
 	result = ((double) totalPairs / (double) mostFreqNr) * ((columnArr.length) * (columnArr.length -1 ) / 2); 
 
 	this.jores = result;
 	
 	}
    
    /**
     * Gets Jores score.
     * 
     * @return Jores score 
     */
    
    double getJores() {
    	
    	if (this.joresInitStat == false) {
    		
    		this.joresScore();
    		
    		return this.jores;
    	}
    	
    	else {
    		
    		return this.jores;
    	}
    }
 

    /**
     * Calculates Schneider score.
     */
	void schneiderScore() {
		
		double normal = 1.0 / 20.0;
		
		double result = ShannonEnthropy.ShannonLn(acidsIntMap, columnArr.length) * normal;
			
		this.schneider = result;
		
		this.schneiderInitStat = true;
	}
	
	 /**
     * Gets Schneider score.
     * 
     * @return Schneider score 
     */
    
    double getSchneider() {
    	
    	if (this.schneiderInitStat == false) {
    		
    		this.schneiderScore();
    		
    		return this.schneider;
    	}
    	
    	else {
    		
    		return this.schneider;
    	}
    }
 
	/**
	 * Calculates Shenkin score.
	 * 
	 * @return Shenkin score
	 */

	void shenkinScore() {
		
		double result = Math.pow( 2.0, ShannonEnthropy.ShannonLog2(acidsIntMap, columnArr.length)) * 6.0;
			
		this.shenkin = result;
		
		this.shenkinInitStat = true;
		
	}
	
	 /**
     * Gets Shenkin score.
     * 
     * @return Shenkin score 
     */
    
    double getShenkin() {
    	
    	if (this.shenkinInitStat == false) {
    		
    		this.shenkinScore();
    		
    		return this.shenkin;
    	}
    	
    	else {
    		
    		return this.shenkin;
    	}
    }
 

	/**
	 * Calculates Gerstein score.
	 * 
	 * @return Gerstein score
	 */
	
	void gersteinScore() {
		
		double result = ShannonEnthropy.ShannonLn(matrix.totalAcidsFrequency(), matrix.numberOfColumns() * matrix.numberOfRows()) - ShannonEnthropy.ShannonLn(acidsIntMap, columnArr.length);
		
		this.gerstein = result;
		
		this.gersteinInitStat = true;
		
	}
	
	 /**
     * Gets Kabat score.
     * 
     * @return Kabat score 
     */
    
    double getGerstein() {
    	
    	if (this.gersteinInitStat == false) {
    		
    		this.gersteinScore();
    		
    		return this.gerstein;
    	}
    	
    	else {
    		
    		return this.gerstein;
    	}
    }
 

 
	/**
	 * Calculates Taylor score. Returns the number of elements in the smallest set covering the whole
	 * column. Gaps sore is included in the largest set. Therefore a column containing gaps is automatically given the lowest possibel score.  
	 *
	 * @return Taylor score
	 */
	
	void SmallestTaylorSetGaps() {
	
	Map<String, HashSet<Character>> setMap = ConservationSets.taylorSets();
	
	Map<String,Integer> repSets = new HashMap<String,Integer>();
	
	Set<String> setMapKeys = setMap.keySet();
	
	Iterator<String> itr = setMapKeys.iterator();
	
	while(itr.hasNext()) {
		
		String key = itr.next();
		
		if (setMap.get(key).containsAll(acidsIntMap.keySet())) {
			
			repSets.put(key, new Integer(setMap.get(key).size()));
				
		}
	}
		
	int smallestSetSize = Collections.min(repSets.values());
	
	assert smallestSetSize > 0 && smallestSetSize < 25;
	
	this.taylorGaps = smallestSetSize;
	
	this.taylorGapsInitStat = true;
	
	}
	
	 /**
     * Gets Kabat score.
     * 
     * @return Kabat score 
     */
    
    double getTaylorGaps() {
    	
    	if (this.taylorGapsInitStat == false) {
    		
    		this.SmallestTaylorSetGaps();
    		
    		return this.taylorGaps;
    	}
    	
    	else {
    		
    		return this.taylorGaps;
    	}
    }
 
 
	/**
	 * Does a very similar thing to SmallestTaylorSetGaps but does not take gaps into account at all.
	 *  
	 * @return Taylor score
	 */
	void SmallestTaylorSetNoGaps() {
	 
	    Map<String, HashSet<Character>> setMap = ConservationSets.taylorSets();
		
		Map<Character,Integer> acidsMapNoGaps = new HashMap<Character,Integer>(acidsIntMap);
		
		if(acidsMapNoGaps.containsKey('-')) {
			
			acidsMapNoGaps.remove('-');
			
		}
		
		Map<String,Integer> repSets = new HashMap<String,Integer>();
		
		Set<String> setMapKeys = setMap.keySet();
		
		Iterator<String> itr = setMapKeys.iterator();
		
		while(itr.hasNext()) {
			
			String key = itr.next();
			
			if (setMap.get(key).containsAll(acidsMapNoGaps.keySet())) {
				
				repSets.put(key, new Integer(setMap.get(key).size()));
					
			}
		}
			
		int smallestSetSize = Collections.min(repSets.values());
		
		assert smallestSetSize > 0 && smallestSetSize < 25;
		
		this.taylorNoGaps = smallestSetSize;
		
		this.taylorNoGapsInitSatat = true;
		
		}
	
	 /**
     * Gets TaylorNoGaps score.
     * 
     * @return TaylorNoGaps score 
     */
    
    double getTaylorNoGaps() {
    	
    	if (this.taylorNoGapsInitSatat == false) {
    		
    		this.SmallestTaylorSetNoGaps();
    		
    		return this.taylorNoGaps;
    	}
    	
    	else {
    		
    		return this.taylorNoGaps;
    	}
    }
 
	/**
	 * Gives a score of one to ten based on whether all the amino acids though the column
	 * maintain or fail to maintain a certain trait.
	 *  
	 * @return
	 */
 	void zvelibilScore(){
 		
 		int result = 0;
 		
 		Map<String, HashSet<Character>> setMap = ConservationSets.zvelibilSets();
 		
 		Set<String> keys = setMap.keySet();
 		
 		Iterator<String> itr = keys.iterator();
 		
 		while(itr.hasNext()) {
 			
 			if(setMap.get(itr.next()).containsAll(acidsIntMap.keySet())) {
 				
 				result++;
 			}
 			
 		}
 		
 		assert result >= 0 && result < 11;
 		
 		this.zvelibil = result;
 		
 		this.zvelibilInitStat = true;
 	}
 	
 	 /**
     * Gets Zvelibil score.
     * 
     * @return Zvelibil score 
     */
    
    double getZvalibil() {
    	
    	if (this.zvelibilInitStat == false) {
    		
    		this.zvelibilScore();
    		
    		return this.zvelibil;
    	}
    	
    	else {
    		
    		return this.zvelibil;
    	}
    }
 

 	/**
 	 * Calculates Karlin score.
 	 *  done exactly like in jon's code: gaps not counted, if all but 1 gaps conservation score 0
 	 *  although it looks to me that the formula suggests that if one doesn't count gaps
 	 *  than the score of all but one gaps would be undefined.
 	 *  
 	 * @return Karlin score
 	 */
 	
 	void karlinScore() {
 		
 		double blosumSum = 0;
 		
 		for( int a = 0; a < columnArr.length; a++) {
 			
 			if(columnArr[a] != '-') {
 			
 				for(int b = a + 1; b < columnArr.length; b++) {
 					
 					if(columnArr[b] != '-') {
 						
 						double pairScore =  ConservationMatrices.BlosumPair(columnArr[a], columnArr[b]);
 						
 						double aSelf =  ConservationMatrices.BlosumPair(columnArr[a], columnArr[a]);
 						
 						assert aSelf > 0;
 						
 						double bSelf = ConservationMatrices.BlosumPair(columnArr[b], columnArr[b]);
 						
 						assert bSelf > 0;
 						
 						blosumSum = blosumSum + ((pairScore) / (Math.sqrt(aSelf * bSelf)));
 					
 					}
 					
 				}
 		
 			}
 			
 		}
 		
 		
 		double finalSum = blosumSum * (2.0 / ( columnArr.length * (columnArr.length - 1)));
 		
 	    assert finalSum >= -1 && finalSum <= 1;
 		
 		this.karlin = finalSum;
 		
 		this.karlinInitStat = true;
 	}
 	
 	 /**
     * Gets Karlin score.
     * 
     * @return Karlin score 
     */
    
    double getKarlin() {
    	
    	if (this.karlinInitStat == false) {
    		
    		this.karlinScore();
    		
    		return this.karlin;
    	}
    	
    	else {
    		
    		return this.karlin;
    	}
    }
 
 	
 	/**
 	 * Calculates Armon score.
 	 * 
 	 * @return Armon score.
 	 */
 // creates an array containing all the amino acids and gaps present in the column
 // iterates through that array twice(nested loops), finds all the possible pairs 
 // that can be formed by aa present
 // gap is considered the 21 aminoacid
 	void armonScore() {
 		
 		double scoreSum = 0;
 		
 		int arrayLength = acidsIntMap.keySet().size();
 		
 		Character[] acidsPresent = new Character[arrayLength];
 		
 		int arrayIndex = 0;
 		
 		Set<Character> keys = acidsIntMap.keySet();
 		
 		Iterator<Character> itr = keys.iterator();
 		
 		while(itr.hasNext()) {
 			
 			acidsPresent[arrayIndex] = itr.next();
 			
 			arrayIndex++;
 			
 		}
 		
 		for ( int a = 0; a < acidsPresent.length; a++ ) {
 			
 			char charA = acidsPresent[a];
 			
 			for (int b = a + 1; b < acidsPresent.length; b++) {
 				
 				char charB = acidsPresent[b];
 				
 				scoreSum = scoreSum + ConservationMatrices.miyataArmonPair(charA, charB);
 				
 			}
 		}
 		
 		this.armon = scoreSum;
 		
 		this.armonInitStat = true;
 	}
 	
 	 /**
     * Gets Armon score.
     * 
     * @return Armon score 
     */
    
    double getArmon() {
    	
    	if (this.armonInitStat == false) {
    		
    		this.armonScore();
    		
    		return this.armon;
    	}
    	
    	else {
    		
    		return this.armon;
    	}
    }
 
 	
 	/** 
 	 * Calculates Thompson score.
 	 * Gaps are accounted for. 
 	 * 
 	 * @return
 	 */
// amino acids are viewed as points in k dimensional space
// an average point is calculated
// score is the distance between the av point and the actual point
 	
 	void thompsonScore(){
 		
 		double[] sum = null;
 		
 		double[] meanPoint = null;
 		
 		double distance = 0.0;
 		
 		double nonGapsFraction = 0.0;
 		
 		char[] alp = Alphabet.alphabetArray();
 		
 		assert alp != null && alp.length != 0;
 		
 		double[][] points = new double[columnArr.length][alp.length];
 		
 		for (int i = 0; i < columnArr.length; i++) {
 			
 			for (int j = 0; j < alp.length; j++) {
 				
 				points[i][j] = ConservationMatrices.BlosumPair(columnArr[i], alp[j]);
 				
 			}
 			
 			if (sum == null) {
 				
 				sum = points[i];
 				
 			}
 				
 				else {
 					
 					sum = ConservationAccessory.addPoints(sum, points[i]);
 				}
 			
 		}
 		
 		assert sum != null;
 		
 		meanPoint = ConservationAccessory.multPointByScalar(sum, 1.0/(double) columnArr.length);
 		
 		for (int i = 0; i < columnArr.length; i++) {
 			
 			distance = distance + ConservationAccessory.pointDistance(points[i], meanPoint);
 			
 		}
 		
 		if (acidsIntMap.keySet().contains('-')) {
 			
 			nonGapsFraction = (double) (columnArr.length - acidsIntMap.get('-')) / (double) columnArr.length;
 			
 		}
 		
 		else {
 			
 			nonGapsFraction = 1.0;
 			
 		}
 		double result = nonGapsFraction * 1.0/columnArr.length * distance;
 		
 		this.thompson  = result;
 		
 		this.thompsonInitStat = true;
 			
 		}
 	
 	 /**
     * Gets Thompson score.
     * 
     * @return Thompson score 
     */
    
    double getThompson() {
    	
    	if (this.thompsonInitStat == false) {
    		
    		this.thompsonScore();
    		
    		return this.thompson;
    	}
    	
    	else {
    		
    		return this.thompson;
    	}
    }
 
 		
 	// causes some math problem because denominator can be 0, that's a formula flaw
 	// nothing can be done about it
 	void lancetScore() {
		
		double result = 0.0;
		
		Set<Character> keys = acidsIntMap.keySet();
		
		Iterator<Character> itr1 = keys.iterator();
		
		while (itr1.hasNext()) {
			
			char key1 = itr1.next();
			
			Iterator<Character> itr2 = keys.iterator();
			
			while (itr2.hasNext()) {
				
				char key2 = itr2.next();
				
				double blosum = ConservationMatrices.BlosumPair(key1, key2);
				
				if (blosum == 0.0) { blosum = 0.00000000000000000000000001;}
				
				result = result + ((((double) acidsIntMap.get(key1)/ (double) columnArr.length * (double) acidsIntMap.get(key1)/ (double) columnArr.length)) / blosum);
			}
		}
	
	this.lancet = result;
	
	this.lancetInitStat  = true;
	
	}
 	
 	 /**
     * Gets Lancet score.
     * 
     * @return Lancet score 
     */
    
    double getLancet() {
    	
    	if (this.lancetInitStat == false) {
    		
    		this.lancetScore();
    		
    		return this.lancet;
    	}
    	
    	else {
    		
    		return this.lancet;
    	}
    }
 
	/**
	 * Calculates Mirny Score.
	 * 
	 * @return mirny score
	 */
 	//calculates shannon enthropy but based on the sets, has to calculate number of amino acids tat belong 
 	//a particular set, stores them in a hashmap
 	//reads in a hashmap with the sets needed and creates a hashmap with set names as keys and number of aa belonging to set as value
 	
	void mirnyScore() {
		
		double mirnySum = 0.0;
		
		Map<String, HashSet<Character>> mirnySets= ConservationSets.mirnySets();
		
		assert mirnySets != null && !mirnySets.isEmpty();
		
		Set<String> mirnyKeys = mirnySets.keySet();
		
		assert !mirnyKeys.isEmpty();
		
		Iterator<String> mirnyKeysItr = mirnyKeys.iterator();
		
		Set<Character> acInKeys = acidsIntMap.keySet(); 
		
		assert !acInKeys.isEmpty();
		
		Map<String,Integer> setsFreq = new HashMap<String,Integer>();  
		
		while (mirnyKeysItr.hasNext()) {
			
			String mirnyKey = mirnyKeysItr.next();
			
			Iterator<Character> acInKeysItr = acInKeys.iterator();
			
			while (acInKeysItr.hasNext()) {
				
				Character acInKey = acInKeysItr.next();
				
				if (mirnySets.get(mirnyKey).contains(acInKey)) {
					
					Integer count = setsFreq.get(mirnyKey);
					
					if (count == null) {
						
						setsFreq.put(mirnyKey, acidsIntMap.get(acInKey));
					}
					
					else {
						
						setsFreq.put(mirnyKey, count + acidsIntMap.get(acInKey));
						
					}
				}
			}
		}
		
		// this assertion will not work if we feed an empty column
		
		assert !setsFreq.isEmpty();
		
		Set<String> setsFreqKeys = setsFreq.keySet();
		
		Iterator<String> setsFreqKeysItr = setsFreqKeys.iterator();
		
		while(setsFreqKeysItr.hasNext()) {
			
			String setFreqKey = setsFreqKeysItr.next();
			
			double pI = (double) setsFreq.get(setFreqKey) / (double) columnArr.length; 
			
			mirnySum = mirnySum + (pI * Math.log(pI));
			
		}
		
		this.mirny = mirnySum;
		
		this.mirnyInitStat = true;
		
	}
	
	 /**
     * Gets Mirny score.
     * 
     * @return Mirny score 
     */
    
    double getMirny() {
    	
    	if (this.mirnyInitStat == false) {
    		
    		this.mirnyScore();
    		
    		return this.mirny;
    	}
    	
    	else {
    		
    		return this.mirny;
    	}
    }
 
	/**
	 * Calculates Williamson score.
	 * 
	 * 
	 * @return Williamson score.
	 */

	void williamsonScore() {
		  
		double willSum = 0.0;
		
		Map<String, HashSet<Character>> willSets= ConservationSets.williamsonSets();
		
		assert willSets != null && ! willSets.isEmpty();
		
		Set<String> willKeys = willSets.keySet();
		
		assert !willKeys.isEmpty();
		
		Iterator<String> willKeysItr = willKeys.iterator();
		
		Set<Character> acInKeys = acidsIntMap.keySet(); 
		
		assert !acInKeys.isEmpty();
		
		Map<String,Integer> setsFreq = new HashMap<String,Integer>();  
		
		while (willKeysItr.hasNext()) {
			
			String willKey = willKeysItr.next();
			
			Iterator<Character> acInKeysItr = acInKeys.iterator();
			
			while (acInKeysItr.hasNext()) {
				
				Character acInKey = acInKeysItr.next();
				
				if (willSets.get(willKey).contains(acInKey)) {
					
					Integer count = setsFreq.get(willKey);
					
					if (count == null) {
						
						setsFreq.put(willKey, acidsIntMap.get(acInKey));
					}
					
					else {
						
						setsFreq.put(willKey, count + acidsIntMap.get(acInKey));
						
					}
				}
			}
		}
		
		// this assertion will not work if we feed an empty column
		
		assert !setsFreq.isEmpty();
		
		Set<String> setsFreqKeys = setsFreq.keySet();
		
		Iterator<String> setsFreqKeysItr = setsFreqKeys.iterator();
		
		while(setsFreqKeysItr.hasNext()) {
			
			String setFreqKey = setsFreqKeysItr.next();
			
			// FIXME Pi in the logarithm needs to be divided by average pi, do it once u get classes to be nested
			
			double pI = (double) setsFreq.get(setFreqKey) / (double) columnArr.length; 
			
			double piAve = (double) matrix.totalAcidsWillSets().get(setFreqKey);
			
			willSum = willSum + (pI * Math.log(pI/piAve));
			
		}
		
		this.williamson = willSum;
		
		this.williamsonInitStat = true;
		
	}
	
	 /**
     * Gets Williamson score.
     * 
     * @return Williamson score 
     */
    
    double getWilliamson() {
    	
    	if (this.williamsonInitStat == false) {
    		
    		this.williamsonScore();
    		
    		return this.williamson;
    	}
    	
    	else {
    		
    		return this.williamson;
    	}
    }
 
	
	/**
	 * Calculates Landgraf score.
	 * 
	 * @return landgraf score
	 */
	void landgrafScore() {
		
		double sum = 0.0;
		
		for (int i = 0; i < columnArr.length; i++) {
			
			for (int j = i + 1; j < columnArr.length; j++) {
				
				double disIJ = ConservationMatrices.dissimilarity(columnArr[i], columnArr[j]);
				
				double disJI = ConservationMatrices.dissimilarity(columnArr[j], columnArr[i]);
				
				sum = sum + (matrix.getVoronoiWeights(1000)[i] * disIJ) + (matrix.getVoronoiWeights(1000)[j] * disJI);
			}
			
		}
		
		double result = sum / columnArr.length;
		
		this.landgraf = result;
		
		 this.landgrafInitStat = true;
	}
	
	 /**
     * Gets Landgraf score.
     * 
     * @return Landgraf score 
     */
    
    double getLandgraf() {
    	
    	if (this.landgrafInitStat == false) {
    		
    		this.landgrafScore();
    		
    		return this.landgraf;
    	}
    	
    	else {
    		
    		return this.landgraf;
    	}
    }
 
	
	/** 
	 * Calculates sander score.
	 * 
	 * @return sander score
	 */ 
	
	void sanderScore() {
		
		double sum = 0.0;
		
		double moderator = 0.0;
		
		for (int i = 0; i < columnArr.length; i++) {
			
			for (int j = i + 1; j < columnArr.length; j++) {
				
				sum = sum + (1 - matrix.getPercentIdentity()[i][j]) * ConservationMatrices.pam250Pair(columnArr[i], columnArr[j]);
			}
			
		}
		
		for (int i = 0; i < columnArr.length; i++) {
			
			for (int j = i + 1; j < columnArr.length; j++) {
			
			moderator = moderator + (1 - matrix.getPercentIdentity()[i][j]);
			
			}
		}

		double result = sum * moderator;
		
		this.sander = result;
		
		this.sanderInitStat = true;
		
		}

	 /**
     * Gets Sander score.
     * 
     * @return Sander score 
     */
    
    double getSander() {
    	
    	if (this.sanderInitStat == false) {
    		
    		this.sanderScore();
    		
    		return this.sander;
    	}
    	
    	else {
    		
    		return this.sander;
    	}
    }
 
	/**
	 * Calculates Valdar score.
	 * 
	 * @return Valdar score
	 */
	void valdarScore() {
		
		double sum = 0.0;
		
		double moderator = 0.0;
		
		for (int i = 0; i < columnArr.length; i++) {
			
			for (int j = i + 1; j < columnArr.length; j++) {
				
				sum = sum + matrix.vingronArgosWeights()[i] * matrix.vingronArgosWeights()[j] * ConservationMatrices.pet91Pair(columnArr[i], columnArr[j]);
			}
			
		}
		
		for (int i = 0; i < columnArr.length; i++) {
			
			for (int j = i + 1; j < columnArr.length; j++) {
			
				moderator = moderator + matrix.vingronArgosWeights()[i] * matrix.vingronArgosWeights()[j];
			}
			
		}
		
		double result = sum * moderator;
		
		this.valdar = result;
		
		this.valdarInitStat = true;
		
		//this.result;
		
	}
	
	 /**
     * Gets Valdar score.
     * 
     * @return Valdar score 
     */
    
    double getValdar() {
    	
    	if (this.valdarInitStat == false) {
    		
    		this.valdarScore();
    		
    		return this.valdar;
    	}
    	
    	else {
    		
    		return this.valdar;
    	}
    }
 
	
} 


