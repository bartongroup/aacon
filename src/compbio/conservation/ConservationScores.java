package compbio.conservation;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

public class ConservationScores {
	
	
	private AminoAcidMatrix matrix;
	
	public ConservationScores(AminoAcidMatrix matrix) {
		
		this.matrix = matrix;
	}
	
	/**
	 * Checks if all but one residues in the column are gaps.
	 * 
	 * @return true if all but one residues are gaps, false if not 
	 */
	
    private boolean allButOneGaps(int columnNr) {
    	
    	if (matrix.getTotalAcidsFreqByCol()[columnNr].containsKey('-') && matrix.getTotalAcidsFreqByCol()[columnNr].get('-') == matrix.getInverseMatrix()[columnNr].length - 1) {
			
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
    private boolean oneResidueTypeNoGaps(int columnNr) {
    	
    	if (matrix.getTotalAcidsFreqByCol()[columnNr].size() == 1 && matrix.getTotalAcidsFreqByCol()[columnNr].containsKey('-') == false) {
    		
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
    
   private boolean containsGaps(int columnNr) {
    	
    	if (matrix.getTotalAcidsFreqByCol()[columnNr].containsKey('-')) {
    		
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
    
   private int numberOfAcidsWithGap(int columnNr) {
    	
    	return matrix.getTotalAcidsFreqByCol()[columnNr].size(); 
    	
    }
    
    /**
     * Counts the number of different amino acids in the column.
     * Gap not counted as a 21 amino acid.
     * 
     * @return the number of different amino acids
     */
    
    private int numberOfAcidsNoGap( int columnNr) {
    
    if(this.containsGaps(columnNr) == true) {
    	
    	return matrix.getTotalAcidsFreqByCol()[columnNr].size() - 1;
    }
    
    else {
    	
    	return matrix.getTotalAcidsFreqByCol()[columnNr].size();
    
    }
    
    }
    /** 
     * Calculates the number of times the most common amino acid in the column occurs.
     * Does not count gap as 21 amino acid.
     * 
     * @return times the most common amino acid occurs
     */
    
    private int mostCommonNumber(int columnNr) {
    	
    	int max = 0;
    	
    	Set<Character> keys = matrix.getTotalAcidsFreqByCol()[columnNr].keySet();
    	
    	Iterator<Character> itr = keys.iterator();
    	
    	while(itr.hasNext()) {
    		
    	Character key = itr.next();
    		
    		if (key != '-' && matrix.getTotalAcidsFreqByCol()[columnNr].get(key) > max) {
    			
    			max = matrix.getTotalAcidsFreqByCol()[columnNr].get(key);
    			
    		}
    	
    	}
    	
    	assert max != 0 : "Zero in the most Common Number";
    	
    	return max;
    	
    }
    
    //Map<Character,Integer> getAcidsIntMap() {
    
    //return acidsIntMap;
    
    
   // }
    /**
     * Gets length of the column.
     * 
     * @return length of the column
     */
    			
    //int length() {
	 
	 //return columnArr.length;
	 
    //}
 
    /**
     * Calculates Kabat score for the column.
     *  
     * @return Kabat score
     */
    double[] kabatScore() {
    	
    	double[] result = new double[matrix.getInverseMatrix().length];
    	
    	for(int i = 0; i < matrix.getInverseMatrix().length; i++) {
	 
        assert mostCommonNumber(i) > 0 && mostCommonNumber(i) < matrix.getInverseMatrix()[i].length + 1;	 
         		
		result[i] = matrix.getInverseMatrix()[i].length * (double) numberOfAcidsNoGap(i)/ (double) mostCommonNumber(i); 
		
    	}
	
		return result;
		
	}
 
    /**
     * Calculates Jores score for the column.
     * When calculating the number of distinct pair formed calculates a combination of two out of the aminoacid types
     * present in the column. Than adds the self pairs. However when there is only one of a type
     * in the column ne self pair is formed. 
     * 
     * @return Jores score 
     */
    double[] joresScore() {
	 
    double[] result	= new double[matrix.getInverseMatrix().length];
    
    Map<Character, Integer> acidsIntMap = null;
    
    Map<Character, Integer> acidsIntMapCopy = null;
    
    Set<Character> keys = null;
    
    Iterator<Integer> itr = null;
    
    Iterator<Character> itr2 = null;
    
    Iterator<Integer> itr3 = null;
    
    	
	for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
	 
	// special case #1 one residue type only
		
	// special case #2 all but one are gaps
	 
	 if (oneResidueTypeNoGaps(i) == true || allButOneGaps(i) == true) {
		 
		 if (oneResidueTypeNoGaps(i) == true) 
				 
				 result[i] = 1.0;
		 
		 else {
			 
			 result[i] = 0.0;
		 }
	 
	 }
	 
	 else {
	 
	 int samePairs = 0;
	 
	 int differentPairs = 0;
	 
	 acidsIntMap = matrix.getTotalAcidsFreqByCol()[i]; 
	 
	 acidsIntMapCopy = new HashMap<Character,Integer>(acidsIntMap);
	 
	 keys = acidsIntMapCopy.keySet();
	 
	 keys.remove('-');
	 
	 int types = keys.size();
	 
	 differentPairs = types * (types - 1) / 2;
	 
	 itr = acidsIntMapCopy.values().iterator();
	 
	 while (itr.hasNext()) {
		 
		 if (itr.next() > 1) {
			 
			 samePairs++;
			 
		 }
		 
	 }
	 
	 int totalPairs = samePairs + differentPairs;
	 
	 itr2 = acidsIntMapCopy.keySet().iterator();
	 
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
	 
	 itr3 = acidsIntMapCopy.values().iterator();
	 
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
		 
	 	}
		 
 	result[i] = ((double) totalPairs / (double) mostFreqNr) * ((matrix.getInverseMatrix()[i].length) * (matrix.getInverseMatrix()[i].length -1 ) / 2); 
 	
	 }
	 
	}
 
 	return result;
 	
 	}
	 
	 
	// Symbol Enthropy Scores

    /**
     * Calculates Schneider score.
     */
    double[] schneiderScore() {
		
		double[] result = new double[matrix.getInverseMatrix().length];
		
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
		
		double normal = 1.0 / Math.log(20.0);
		
		result[i] = ShannonEnthropy.ShannonLn(matrix.getTotalAcidsFreqByCol()[i], matrix.getInverseMatrix()[i].length) * normal;
		
		assert result[i] >= 0 && result[i] <= 1;
		
		}
			
		return result;
	}
	
	/**
	 * Calculates Shenkin score.
	 * 
	 * @return Shenkin score
	 */

	double[] shenkinScore() {
		
		double[] result = new double[matrix.getInverseMatrix().length];
		
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
		
		result[i] = Math.pow( 2.0, ShannonEnthropy.ShannonLog2(matrix.getTotalAcidsFreqByCol()[i], matrix.getInverseMatrix()[i].length)) * 6.0;
			
		assert result[i] >= 6 && result[i] <= 120;
		
		}
		
		return result;
	}

	/**
	 * Calculates Gerstein score.
	 * 
	 * @return Gerstein score
	 */
	
	double[] gersteinScore() {
		
		double[] result = new double[matrix.getInverseMatrix().length];
		
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
		
		result[i] = ShannonEnthropy.ShannonLn(matrix.totalAcidsFrequency(), matrix.numberOfColumns() * matrix.numberOfRows()) - ShannonEnthropy.ShannonLn(matrix.getTotalAcidsFreqByCol()[i], matrix.getInverseMatrix()[i].length);
		
		}
		
		return result;
		
	}

 
	/**
	 * Calculates Taylor score. Returns the number of elements in the smallest set covering the whole
	 * column. Gaps sore is included in the largest set. Therefore a column containing gaps is automatically given the lowest possibel score.  
	 *
	 * @return Taylor score
	 */
	
	double[] SmallestTaylorSetGaps() {
	
	Map<String, HashSet<Character>> setMap = ConservationSets.taylorSets();
	
	double[] smallestSetSize = new double[matrix.getInverseMatrix().length];
	
	Map<String,Integer> repSets = null;
	
	Set<String> setMapKeys = null;
	
	Iterator<String> itr = null;
	
	for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
	
	repSets = new HashMap<String,Integer>();
	
	setMapKeys = setMap.keySet();
	
	itr = setMapKeys.iterator();
	
	while(itr.hasNext()) {
		
		String key = itr.next();
		
		if (setMap.get(key).containsAll(matrix.getTotalAcidsFreqByCol()[i].keySet())) {
			
			repSets.put(key, new Integer(setMap.get(key).size()));
				
		}
	}
		
	smallestSetSize[i] = Collections.min(repSets.values());
	
	assert smallestSetSize[i] > 0 && smallestSetSize[i] < 25;
	
	}
	
	return smallestSetSize;
	
	}
 
	/**
	 * Does a very similar thing to SmallestTaylorSetGaps but does not take gaps into account at all.
	 *  
	 * @return Taylor score
	 */
	double[] SmallestTaylorSetNoGaps() {
	 
	    Map<String, HashSet<Character>> setMap = ConservationSets.taylorSets();
	    
	    Set<String> setMapKeys = setMap.keySet();
	    
	    double[] smallestSetSize = new double[matrix.getInverseMatrix().length];
	    
	    Map<Character,Integer> acidsMapNoGaps = null;
	    
	    Map<String,Integer> repSets = null;
		
		Iterator<String> itr = null;
	    
	    for(int i = 0; i < matrix.getInverseMatrix().length; i++) {
		
		acidsMapNoGaps = new HashMap<Character,Integer>(matrix.getTotalAcidsFreqByCol()[i]);
		
		if(acidsMapNoGaps.containsKey('-')) {
			
			acidsMapNoGaps.remove('-');
			
		}
		
		repSets = new HashMap<String,Integer>();
		
		itr = setMapKeys.iterator();
		
		while(itr.hasNext()) {
			
			String key = itr.next();
			
			if (setMap.get(key).containsAll(acidsMapNoGaps.keySet())) {
				
				repSets.put(key, new Integer(setMap.get(key).size()));
					
			}
		}
			
		smallestSetSize[i] = Collections.min(repSets.values());
		
		assert smallestSetSize[i] > 0 && smallestSetSize[i] < 25;
		
	    }
		
		return smallestSetSize;
		
		}
	/**
	 * Gives a score of one to ten based on whether all the amino acids though the column
	 * maintain or fail to maintain a certain trait.
	 *  
	 * @return
	 */
 	double[] zvelibilScore(){
 		
 		double[] result = new double[matrix.getInverseMatrix().length];
 		
 		Map<String, HashSet<Character>> setMap = ConservationSets.zvelibilSets();
 		
 		Set<String> keys = setMap.keySet();
 		
 		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
 		
 		Iterator<String> itr = keys.iterator();
 		
 		while(itr.hasNext()) {
 			
 			if(setMap.get(itr.next()).containsAll(matrix.getTotalAcidsFreqByCol()[i].keySet())) {
 				
 				result[i]++;
 			}
 			
 		}
 		
 		assert result[i] >= 0 && result[i] < 11;
 		
 		}
 		
 		return result;
 	}

 	/**
 	 * Calculates Karlin score.
 	 *  done exactly like in jon's code: gaps not counted, if all but 1 gaps conservation score 0
 	 *  although it looks to me that the formula suggests that if one doesn't count gaps
 	 *  than the score of all but one gaps would be undefined.
 	 *  
 	 * @return Karlin score
 	 */
 	
 	double[] karlinScore() {
 		
 		double[] finalSum = new double[matrix.getInverseMatrix().length];
 		
 		double[] blosumSum = new double[matrix.getInverseMatrix().length];
 		
 		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
 		
 		for( int a = 0; a < matrix.getInverseMatrix()[i].length; a++) {
 			
 			if(matrix.getInverseMatrix()[i][a] != '-') {
 			
 				for(int b = a + 1; b < matrix.getInverseMatrix()[i].length; b++) {
 					
 					if(matrix.getInverseMatrix()[i][b] != '-') {
 						
 						double pairScore =  ConservationMatrices.BlosumPair(matrix.getInverseMatrix()[i][a], matrix.getInverseMatrix()[i][b]);
 						
 						double aSelf =  ConservationMatrices.BlosumPair(matrix.getInverseMatrix()[i][a], matrix.getInverseMatrix()[i][a]);
 						
 						assert aSelf > 0;
 						
 						double bSelf = ConservationMatrices.BlosumPair(matrix.getInverseMatrix()[i][b], matrix.getInverseMatrix()[i][b]);
 						
 						assert bSelf > 0;
 						
 						blosumSum[i] = blosumSum[i] + ((pairScore) / (Math.sqrt(aSelf * bSelf)));
 					
 					}
 					
 				}
 		
 			}
 			
 		}
 		
 		finalSum[i] = blosumSum[i] * (2.0 / ( matrix.getInverseMatrix()[i].length * (matrix.getInverseMatrix()[i].length - 1)));
 		
 	    assert finalSum[i] >= -1 && finalSum[i] <= 1;
 	    
 	}
 		
 		return finalSum;
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
 	double[] armonScore() {
 		
 		double[] scoreSum = new double[matrix.getInverseMatrix().length];
 		
 		Set<Character> keys = null;
 		
 		Iterator<Character> itr = null;
 		
 		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
 		
 		int arrayLength = matrix.getTotalAcidsFreqByCol()[i].keySet().size();
 		
 		Character[] acidsPresent = new Character[arrayLength];
 		
 		int arrayIndex = 0;
 		
 		keys = matrix.getTotalAcidsFreqByCol()[i].keySet();
 		
 		itr = keys.iterator();
 		
 		while(itr.hasNext()) {
 			
 			acidsPresent[arrayIndex] = itr.next();
 			
 			arrayIndex++;
 			
 		}
 		
 		for ( int a = 0; a < acidsPresent.length; a++ ) {
 			
 			char charA = acidsPresent[a];
 			
 			for (int b = a + 1; b < acidsPresent.length; b++) {
 				
 				char charB = acidsPresent[b];
 				
 				scoreSum[i] = scoreSum[i] + ConservationMatrices.miyataArmonPair(charA, charB);
 				
 			}
 		}
 		
 		}
 		
 		return scoreSum;
 		
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
 	
 	double[] thompsonScore(){
 		
 		double[] sum = null;
 		
 		double[] meanPoint = null;
 		
 		double distance = 0.0;
 		
 		double nonGapsFraction = 0.0;
 		
 		char[] alp = Alphabet.alphabetArray();
 		
 		assert alp != null && alp.length != 0;
 		
 		double[] result = new double[matrix.getInverseMatrix().length]; 
 		
 		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
 		
 		double[][] points = new double[matrix.getInverseMatrix()[i].length][alp.length];
 		
 		for (int a = 0; a < matrix.getInverseMatrix()[i].length; a++) {
 			
 			for (int b = 0; b < alp.length; b++) {
 				
 				points[a][b] = ConservationMatrices.BlosumPair(matrix.getInverseMatrix()[i][a], alp[b]);
 				
 			}
 			
 			if (sum == null) {
 				
 				sum = points[a];
 				
 			}
 				
 				else {
 					
 					sum = ConservationAccessory.addPoints(sum, points[a]);
 				}
 			
 		}
 		
 		assert sum != null;
 		
 		meanPoint = ConservationAccessory.multPointByScalar(sum, 1.0/(double) matrix.getInverseMatrix()[i].length);
 		
 		for (int c = 0; c < matrix.getInverseMatrix()[i].length; c++) {
 			
 			distance = distance + ConservationAccessory.pointDistance(points[c], meanPoint);
 			
 		}
 		
 		if (matrix.getTotalAcidsFreqByCol()[i].keySet().contains('-')) {
 			
 			nonGapsFraction = (double) (matrix.getInverseMatrix()[i].length - matrix.getTotalAcidsFreqByCol()[i].get('-')) / (double) matrix.getInverseMatrix()[i].length;
 			
 		}
 		
 		else {
 			
 			nonGapsFraction = 1.0;
 			
 		}
 		
 		result[i] = nonGapsFraction * 1.0/matrix.getInverseMatrix()[i].length * distance;
 		
 		sum = null;
 		
 		meanPoint = null;
 		
 		distance = 0.0;
 		
 		nonGapsFraction = 0.0;
 		
 		}
 		
 		return result;
 			
 		}
 		
 	// causes some math problem because denominator can be 0, that's a formula flaw
 	// nothing can be done about it
 	double[] notLancetScore() {
		
		double[] result = new double[matrix.getInverseMatrix().length];
		
		Set<Character> keys = null;
		
		Iterator<Character> itr1 = null;
		
		Iterator<Character> itr2 = null;
		
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
		
		keys = matrix.getTotalAcidsFreqByCol()[i].keySet();
		
		itr1 = keys.iterator();
		
		while (itr1.hasNext()) {
			
			char key1 = itr1.next();
			
			itr2 = keys.iterator();
			
			while (itr2.hasNext()) {
				
				char key2 = itr2.next();
				
				double blosum = ConservationMatrices.BlosumPair(key1, key2);
				
				result[i] = result[i] + ((((double)matrix.getTotalAcidsFreqByCol()[i].get(key1)/ (double) matrix.getInverseMatrix()[i].length * (double) matrix.getTotalAcidsFreqByCol()[i].get(key2)/ (double) matrix.getInverseMatrix()[i].length)) * blosum);
			}
		}
		
		}
	
	return result;
	
	}
 	
	/**
	 * Calculates Mirny Score.
	 * 
	 * @return mirny score
	 */
 	//calculates shannon enthropy but based on the sets, has to calculate number of amino acids tat belong 
 	//a particular set, stores them in a hashmap
 	//reads in a hashmap with the sets needed and creates a hashmap with set names as keys and number of aa belonging to set as value
 	
	double[] mirnyScore() {
		
		double[] mirnySum = new double[matrix.getInverseMatrix().length];
		
		Map<String, HashSet<Character>> mirnySets= ConservationSets.mirnySets();
		
		assert mirnySets != null && !mirnySets.isEmpty();
		
		Set<String> mirnyKeys = mirnySets.keySet();
		
		assert !mirnyKeys.isEmpty();
		
		Iterator<String> mirnyKeysItr = null;
		
		Set<Character> acInKeys = null;
		
		Map<String,Integer> setsFreq = null;
		
		Iterator<Character> acInKeysItr = null;
		
		Set<String> setsFreqKeys = null;
		
		Iterator<String> setsFreqKeysItr = null;
		
		for (int i = 0; i < matrix.getInverseMatrix().length; i++ ) {
		
		mirnyKeysItr = mirnyKeys.iterator();
		
		acInKeys = matrix.getTotalAcidsFreqByCol()[i].keySet(); 
		
		assert !acInKeys.isEmpty();
		
		setsFreq = new HashMap<String,Integer>();  
		
		while (mirnyKeysItr.hasNext()) {
			
			String mirnyKey = mirnyKeysItr.next();
			
			acInKeysItr = acInKeys.iterator();
			
			while (acInKeysItr.hasNext()) {
				
				Character acInKey = acInKeysItr.next();
				
				if (mirnySets.get(mirnyKey).contains(acInKey)) {
					
					Integer count = setsFreq.get(mirnyKey);
					
					if (count == null) {
						
						setsFreq.put(mirnyKey, matrix.getTotalAcidsFreqByCol()[i].get(acInKey));
					}
					
					else {
						
						setsFreq.put(mirnyKey, count + matrix.getTotalAcidsFreqByCol()[i].get(acInKey));
						
					}
				}
			}
		}
		
		// this assertion will not work if we feed an empty column
		
		assert !setsFreq.isEmpty();
		
		setsFreqKeys = setsFreq.keySet();
		
		setsFreqKeysItr = setsFreqKeys.iterator();
		
		while(setsFreqKeysItr.hasNext()) {
			
			String setFreqKey = setsFreqKeysItr.next();
			
			double pI = (double) setsFreq.get(setFreqKey) / (double) matrix.getInverseMatrix()[i].length; 
			
			mirnySum[i] = mirnySum[i] + (pI * Math.log(pI));
			
		}
		
		}
		
		return mirnySum;
		
	}
	/**
	 * Calculates Williamson score.
	 * 
	 * 
	 * @return Williamson score.
	 */

	double[] williamsonScore() {
		  
		double[] willSum = new double[matrix.getInverseMatrix().length];
		
		Map<String, HashSet<Character>> willSets= ConservationSets.williamsonSets();
		
		assert willSets != null && ! willSets.isEmpty();
		
		Set<String> willKeys = willSets.keySet();
		
		assert !willKeys.isEmpty();
		
		Iterator<String> willKeysItr = null;
		
		Set<Character> acInKeys = null;
		
		Map<String,Integer> setsFreq = null;
		
		Iterator<Character> acInKeysItr = null;
		
		Set<String> setsFreqKeys = null;
		
		Iterator<String> setsFreqKeysItr = null;
		
		
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
		
		willKeysItr = willKeys.iterator();
		
		acInKeys = matrix.getTotalAcidsFreqByCol()[i].keySet(); 
		
		assert !acInKeys.isEmpty();
		
		setsFreq = new HashMap<String,Integer>();  
		
		while (willKeysItr.hasNext()) {
			
			String willKey = willKeysItr.next();
			
			acInKeysItr = acInKeys.iterator();
			
			while (acInKeysItr.hasNext()) {
				
				Character acInKey = acInKeysItr.next();
				
				if (willSets.get(willKey).contains(acInKey)) {
					
					Integer count = setsFreq.get(willKey);
					
					if (count == null) {
						
						setsFreq.put(willKey, matrix.getTotalAcidsFreqByCol()[i].get(acInKey));
					}
					
					else {
						
						setsFreq.put(willKey, count + matrix.getTotalAcidsFreqByCol()[i].get(acInKey));
						
					}
				}
			}
		}
		
		// this assertion will not work if we feed an empty column
		
		assert !setsFreq.isEmpty();
		
		setsFreqKeys = setsFreq.keySet();
		
		setsFreqKeysItr = setsFreqKeys.iterator();
		
		while(setsFreqKeysItr.hasNext()) {
			
			String setFreqKey = setsFreqKeysItr.next();
			
			// FIXME Pi in the logarithm needs to be divided by average pi, do it once u get classes to be nested
			
			double pI = (double) setsFreq.get(setFreqKey) / (double) matrix.getInverseMatrix()[i].length; 
			
			double piAve = (double) matrix.totalAcidsWillSets().get(setFreqKey);
			
			willSum[i] = willSum[i] + (pI * Math.log(pI/piAve));
			
		}
		
		}
		
		return willSum;
		
	}
	
	/**
	 * Calculates Landgraf score.
	 * 
	 * @return landgraf score
	 */
	double[] landgrafScore() {
		
		double sum = 0;
		
		double[] result = new double[matrix.getInverseMatrix().length];
		
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
		
		for (int a = 0; a < matrix.getInverseMatrix()[i].length; a++) {
			
			for (int b = a + 1; b < matrix.getInverseMatrix()[i].length; b++) {
				
				double disIJ = ConservationMatrices.dissimilarity(matrix.getInverseMatrix()[i][a], matrix.getInverseMatrix()[i][b]);
				
				double disJI = ConservationMatrices.dissimilarity(matrix.getInverseMatrix()[i][b], matrix.getInverseMatrix()[i][a]);
				
				sum = sum + (matrix.getVoronoiWeights(1000)[a] * disIJ) + (matrix.getVoronoiWeights(1000)[b] * disJI);
			}
			
		}
		
		result[i] = sum / matrix.getInverseMatrix()[i].length;
		
		sum = 0;
		
		}
		
		return result;
	}
	
	/** 
	 * Calculates sander score.
	 * 
	 * @return sander score
	 */ 
	
	double[] sanderScore() {
		
		double sum = 0.0;
		
		double moderator = 0.0;
		
		double[] result = new double[matrix.getInverseMatrix().length]; 
		
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			
		for (int a = 0; a < matrix.getInverseMatrix()[i].length; a++) {
			
			for (int b = a + 1; b < matrix.getInverseMatrix()[i].length; b++) {
				
				sum = sum + (1 - matrix.getPercentIdentity()[a][b]) * ConservationMatrices.pam250Pair(matrix.getInverseMatrix()[i][a], matrix.getInverseMatrix()[i][b]);
			}
			
		}
		
		for (int a = 0; a < matrix.getInverseMatrix()[i].length; a++) {
			
			for (int b = a + 1; b < matrix.getInverseMatrix()[i].length; b++) {
			
			moderator = moderator + (1 - matrix.getPercentIdentity()[a][b]);
			
			}
		}

		result[i] = sum * moderator;
		
		sum = 0.0;
		
		moderator = 0.0;
		
		}
		
		return result;
		
		}

	/**
	 * Calculates Valdar score.
	 * 
	 * @return Valdar score
	 */
	double[] valdarScore() {
		
		double sum = 0.0;
		
		double moderator = 0.0;
		
		double[] result = new double[matrix.getInverseMatrix().length];
		
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
		
		for (int a = 0; a < matrix.getInverseMatrix()[i].length; a++) {
			
			for (int b = a + 1; b < matrix.getInverseMatrix()[i].length; b++) {
				
				sum = sum + matrix.vingronArgosWeights()[a] * matrix.vingronArgosWeights()[b] * ConservationMatrices.pet91Pair(matrix.getInverseMatrix()[i][a], matrix.getInverseMatrix()[i][b]);
			}
			
		}
		
		for (int a = 0; a < matrix.getInverseMatrix()[i].length; a++) {
			
			for (int b = a + 1; b < matrix.getInverseMatrix()[i].length; b++) {
			
				moderator = moderator + matrix.vingronArgosWeights()[a] * matrix.vingronArgosWeights()[b];
			}
			
		}
		
		result[i] = sum * moderator;
		
		sum = 0.0;
		
		moderator= 0.0;
		
		}
		
		return result;
		
	}
	

}
