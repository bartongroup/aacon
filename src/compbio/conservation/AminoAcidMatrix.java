package compbio.conservation;

//creates a matrix of aa in multiple alignment
//gets fasta sequences and puts them into matrix 
//might have to check if all sequences are equal length(ask)

import java.util.*;

import compbio.util.SequenceUtil;
import compbio.util.FastaSequence;
import java.io.InputStream;
import java.io.IOException;

/** 
 *  This class provides representation of an alignment as a matrix implemented as 2D array.
 *  Rows correspond to the sequences. Columns correspond to the vertical columns in the alignment consisting of amino acids with the same index in all the   sequences.
 *  The only condition is that all the the sequences in the fasta file are  of the   same length.
 * 	
 * @author agolicz
 *
 */


 class AminoAcidMatrix {
	
	/** 
	 * Stores the matrix.
	 */
	
	private final char[][] matrix;
	
	/**
	 * Stores inverse matrix
	 */
	
	private final char[][] inverseMatrix;
	
	/**
	 * Stores occurances of amino acids in columns, columns indexed starting fromm 0
	 */
	
	private Map<Character, Integer>[] acidsIntMap; 
	
	/**
	 * Holds the in the indices of Xs changed to gaps.
	 * The row number is a key and the column number is value.
	 */
	
	private List<HashMap<Integer, Integer>> xToGapSubs = null;
	
	/**
	 * The total occurrence of amino acids in the whole alignment.
	 */
	
	private Map<Character,Integer> totalFrequency = null;
	
	/**
	 * The total number of amino acids in the whole alignment belonging to each Williamson Set.
	 */
	
	private Map<String, Integer> willSetsTotal = null;
	
	/** 
	 * Vingron Argos weights of the the sequences. 
	 */
	
	private double[] vingronArgosWeights = null;
	
	/**
	 * Percent identity.
	 */
	
	private double[][] percentIdentity = null;
	
	/**
	 * Weights according to Voronoi.
	 */
	private double[] voronoiWeighths = null;
	
	/** 
	 * This constructor constructor allows manual creation of only one column.
	 * Might be of help if somebody wants to check the the functionality of the class without feeding it the whole alignment.
	 * 
	 * @param p1
	 * @param p2
	 * @param p3
	 * @param p4
	 * @param p5
	 * @param p6
	 * @param p7
	 * @param p8
	 * @param p9
	 * @param p10
	 */
	
	public AminoAcidMatrix(char p1, char p2, char p3, char p4, char p5, char p6,char p7 , char p8 ,char p9 ,char p10){
		
		matrix = new char[10][1];
		                        
		matrix[0][0] = p1;
		matrix[1][0] = p2;
		matrix[2][0] = p3;
		matrix[3][0] = p4;
		matrix[4][0] = p5;
		matrix[5][0] = p6;
		matrix[6][0] = p7;
		matrix[7][0] = p8;
		matrix[8][0] = p9;
		matrix[9][0] = p10;
		
		inverseMatrix = null;
		
	}
	
	/**
	 * This constructor allows manual creation of a 10 by 10 matrix.
	 * @param p1
	 * @param p2
	 * @param p3
	 * @param p4
	 * @param p5
	 * @param p6
	 * @param p7
	 * @param p8
	 * @param p9
	 * @param p10
	 * @param p11
	 * @param p12
	 * @param p13
	 * @param p14
	 * @param p15
	 * @param p16
	 * @param p17
	 * @param p18
	 * @param p19
	 * @param p20
	 * @param p21
	 * @param p22
	 * @param p23
	 * @param p24
	 * @param p25
	 * @param p26
	 * @param p27
	 * @param p28
	 * @param p29
	 * @param p30
	 * @param p31
	 * @param p32
	 * @param p33
	 * @param p34
	 * @param p35
	 * @param p36
	 * @param p37
	 * @param p38
	 * @param p39
	 * @param p40
	 * @param p41
	 * @param p42
	 * @param p43
	 * @param p44
	 * @param p45
	 * @param p46
	 * @param p47
	 * @param p48
	 * @param p49
	 * @param p50
	 * @param p51
	 * @param p52
	 * @param p53
	 * @param p54
	 * @param p55
	 * @param p56
	 * @param p57
	 * @param p58
	 * @param p59
	 * @param p60
	 * @param p61
	 * @param p62
	 * @param p63
	 * @param p64
	 * @param p65
	 * @param p66
	 * @param p67
	 * @param p68
	 * @param p69
	 * @param p70
	 * @param p71
	 * @param p72
	 * @param p73
	 * @param p74
	 * @param p75
	 * @param p76
	 * @param p77
	 * @param p78
	 * @param p79
	 * @param p80
	 * @param p81
	 * @param p82
	 * @param p83
	 * @param p84
	 * @param p85
	 * @param p86
	 * @param p87
	 * @param p88
	 * @param p89
	 * @param p90
	 * @param p91
	 * @param p92
	 * @param p93
	 * @param p94
	 * @param p95
	 * @param p96
	 * @param p97
	 * @param p98
	 * @param p99
	 * @param p100
	 */
	
	public AminoAcidMatrix(char p1, char p2, char p3, char p4, char p5, char p6,char p7 , char p8 ,char p9 ,char p10, char p11, char p12, char p13, char p14, char p15, char p16, char p17 ,char p18 ,char p19 ,char p20,char p21, char p22, char p23, char p24, char p25, char p26, char p27 ,char p28 ,char p29 ,char p30,char p31, char p32, char p33, char p34, char p35, char p36,char p37 , char p38 ,char p39 ,char p40,char p41, char p42, char p43, char p44, char p45, char p46,char p47 , char p48 ,char p49 ,char p50,char p51, char p52, char p53, char p54, char p55, char p56,char p57 , char p58 ,char p59 ,char p60, char p61, char p62, char p63, char p64, char p65, char p66,char p67 , char p68 ,char p69 ,char p70,char p71, char p72, char p73, char p74, char p75, char p76,char p77 , char p78 ,char p79 ,char p80,char p81, char p82, char p83, char p84, char p85, char p86,char p87 , char p88 ,char p89 ,char p90,char p91, char p92, char p93, char p94, char p95, char p96,char p97 , char p98 ,char p99 ,char p100){
		
		matrix = new char[10][10];
		                        
		matrix[0][0] = p1;
		matrix[0][1] = p2;
		matrix[0][2] = p3;
		matrix[0][3] = p4;
		matrix[0][4] = p5;
		matrix[0][5] = p6;
		matrix[0][6] = p7;
		matrix[0][7] = p8;
		matrix[0][8] = p9;
		matrix[0][9] = p10;
		matrix[1][0] = p11;
		matrix[1][1] = p12;
		matrix[1][2] = p13;
		matrix[1][3] = p14;
		matrix[1][4] = p15;
		matrix[1][5] = p16;
		matrix[1][6] = p17;
		matrix[1][7] = p18;
		matrix[1][8] = p19;
		matrix[1][9] = p20;
		matrix[2][0] = p21;
		matrix[2][1] = p22;
		matrix[2][2] = p23;
		matrix[2][3] = p24;
		matrix[2][4] = p25;
		matrix[2][5] = p26;
		matrix[2][6] = p27;
		matrix[2][7] = p28;
		matrix[2][8] = p29;
		matrix[2][9] = p30;
		matrix[3][0] = p31;
		matrix[3][1] = p32;
		matrix[3][2] = p33;
		matrix[3][3] = p34;
		matrix[3][4] = p35;
		matrix[3][5] = p36;
		matrix[3][6] = p37;
		matrix[3][7] = p38;
		matrix[3][8] = p39;
		matrix[3][9] = p40;
		matrix[4][0] = p41;
		matrix[4][1] = p42;
		matrix[4][2] = p43;
		matrix[4][3] = p44;
		matrix[4][4] = p45;
		matrix[4][5] = p46;
		matrix[4][6] = p47;
		matrix[4][7] = p48;
		matrix[4][8] = p49;
		matrix[4][9] = p50;
		matrix[5][0] = p51;
		matrix[5][1] = p52;
		matrix[5][2] = p53;
		matrix[5][3] = p54;
		matrix[5][4] = p55;
		matrix[5][5] = p56;
		matrix[5][6] = p57;
		matrix[5][7] = p58;
		matrix[5][8] = p59;
		matrix[5][9] = p60;
		matrix[6][0] = p61;
		matrix[6][1] = p62;
		matrix[6][2] = p63;
		matrix[6][3] = p64;
		matrix[6][4] = p65;
		matrix[6][5] = p66;
		matrix[6][6] = p67;
		matrix[6][7] = p68;
		matrix[6][8] = p69;
		matrix[6][9] = p70;
		matrix[7][0] = p71;
		matrix[7][1] = p72;
		matrix[7][2] = p73;
		matrix[7][3] = p74;
		matrix[7][4] = p75;
		matrix[7][5] = p76;
		matrix[7][6] = p77;
		matrix[7][7] = p78;
		matrix[7][8] = p79;
		matrix[7][9] = p80;
		matrix[8][0] = p81;
		matrix[8][1] = p82;
		matrix[8][2] = p83;
		matrix[8][3] = p84;
		matrix[8][4] = p85;
		matrix[8][5] = p86;
		matrix[8][6] = p87;
		matrix[8][7] = p88;
		matrix[8][8] = p89;
		matrix[8][9] = p90;
		matrix[9][0] = p91;
		matrix[9][1] = p92;
		matrix[9][2] = p93;
		matrix[9][3] = p94;
		matrix[9][4] = p95;
		matrix[9][5] = p96;
		matrix[9][6] = p97;
		matrix[9][7] = p98;
		matrix[9][8] = p99;
		matrix[9][9] = p100;
		
		inverseMatrix = null;
	}
	
	/**
	 * Constructor that reads in a fasta file and creates an amino acid matrix.
	 * Gaps indicated by any sign are now indicated by a dash.
	 * The unknown amino acid X is replaced by a dash.
	 * 
	 * @param inStream
	 */
		
	public AminoAcidMatrix(List<FastaSequence> seqs){
		
			  Set<Character> alph = Alphabet.alphabet();
	          
	          int sequenceNr = seqs.size();
	           
	          FastaSequence seq = seqs.get(0);
	     
	          String firstSequence = seq.getSequence();

	          int sequenceLength = firstSequence.length();
	       
	          matrix = new char[sequenceNr][sequenceLength];
	          
	          inverseMatrix = new char[sequenceLength][sequenceNr];

	          for( int i = 0; i < sequenceNr; i++) {

	                FastaSequence s = seqs.get(i);
	         
	                char[] sequenceChars = s.getSequence().toCharArray();
	                   
			              for ( int j = 0; j < sequenceLength; j++) {
			            	  
			            	  		char ch = sequenceChars[j];
			            	  		
			            	  		if(ch == '.' || ch == '*' || ch == ' ' || ch =='X') {
			        		    		
			        					ch = '-';
			        				}
			            	  		
			            	  		if (ch == 'X') {
			            	  			
			            	  			this.xToGapSubs.add(new HashMap<Integer, Integer>());
			            	  			
			            	  			this.xToGapSubs.get(this.xToGapSubs.size() - 1).put(i, j);
			            	  			
			            	  			ch = '-';
			            	  			
			            	  		}
			        	
			        				assert alph.contains(ch) : "Illegal character in the matrix";
			        				
				
	                                 matrix[i][j] = ch;
	                                 
	                                 inverseMatrix[j][i] = ch;
	                        
	                      }
	          }
	          
	          this.calTotalAcidsFreqByCol();
	          
	}           

	/**
	 * Gets the number of columns.
	 * 
	 * @return number of columns
	 * 
	 */
	
	int numberOfColumns() {

	int nrColumns = matrix[0].length;

	return nrColumns;

	}
	/**
	 * Gets the number of rows.
	 * 
	 * @return number of rows
	 */
	
	int numberOfRows() {

	int nrRows = matrix.length;

	return nrRows;

	}

	/** 
	 * Gets the whole matrix
	 * 
	 * @return matrix
	 */
	char[][] getMatrix() {

	return matrix;

	}
	
	/**
	 * Gets the inverse matrix
	 * 
	 * @return inverseMatrix
	 */
	
	char[][] getInverseMatrix() {
		
		return inverseMatrix;
	}

	/** 
	 * Returns character at a specified position.
	 * 
	 * @param row
	 * @param column
	 * @return character
	 */
	
	char getMatrixPosition(int row, int column) {

	char position = matrix[row][column];

	return position;

	}
	
	/** 
	 * Gets a column specified by given number. Indexing starts with 0.
	 * 
	 * @param number
	 * @return column as an array of characters
	 */
	
	char[] getColumn(int number) {
		
		assert number < this.numberOfColumns();

		Set<Character> alph = Alphabet.alphabet();
		
		char[] column = new char[this.numberOfRows()];

		for (int i = 0; i < this.numberOfRows(); i++) {

			char ch = matrix[i][number];

			if(ch == '.' || ch == '*' || ch == ' ' || ch =='X') {
    		
    			ch = '-';

    			}
    	
    			assert alph.contains(ch) : "Illegal character in the column";

			column[i] = ch;
		}
		
		return column;
		
	}
	
	/**
	 * Gets row specified by a given number. Indexing starts with 0.
	 * @param number
	 * @return row as an array of characters
	 */
	
	char[] getRow(int number) {
		
		assert number < this.numberOfRows();
		
		char[] row = matrix[number];
		
		return row;
		
	}
	
	private void calTotalAcidsFreqByCol() {
		
		acidsIntMap = new Map[this.inverseMatrix.length];
		
		for (int i = 0; i < this.numberOfColumns(); i++) {
			
			acidsIntMap[i] = Alphabet.calculateOccurance(this.inverseMatrix[i]);
		    
		}
		
	}
	
	Map<Character, Integer>[] getTotalAcidsFreqByCol() {
		
		if(acidsIntMap == null) {
			
			this.calTotalAcidsFreqByCol();
		}
		
		return acidsIntMap;
	}
	
	/**
	 * Gets the total amino acid occurrence in the whole alignment.
	 *  
	 * @return map with characters as keys and occurrence as values.
	 */
	Map<Character, Integer> totalAcidsFrequency() {
		
		if (totalFrequency == null) {
			
			this.calTotalAcidsFrequency();
		}
		
		return totalFrequency;
	}
	
	/** 
	 * Gets the total number of amino acids belonging to particular Williamson sets.
	 * 
	 * @return map with set names as keys and number of amino acids belonging to the set as value.
	 * 
	 */
	
	Map<String,Integer> totalAcidsWillSets() {
		
		if(willSetsTotal == null) {
			
			this.calTotalAcidsWillSets();
		}
	
		return willSetsTotal;
	}
	
	/**
	 * Calculates the total amino acid occurrence in the whole alignment.
	 */
	
	    private void calTotalAcidsFrequency() {
		
		Map<Character,Integer> totalFreq = new HashMap<Character,Integer>();
		
		Set<Character> alph = Alphabet.alphabet();
		
		for (int i = 0; i < this.numberOfRows(); i++) {
			
			for (int j = 0; j < this.numberOfColumns(); j++) {
				
				Character ch = matrix[i][j];
				
				if(ch == '.' || ch == '*' || ch == ' ' || ch =='X') {
		    		
					ch = '-';
				}
	
				assert alph.contains(ch) : "Illegal character in the matrix";
				
				Integer count = totalFreq.get(ch);
				
				if (count == null) {
					
					totalFreq.put(ch, 1);
				}
				
				else {
					
					totalFreq.put(ch, count + 1);
				}
			}
		}
		
		totalFrequency = totalFreq;
	}
	    
	    /**
	     * Calculates the total number of amino acids belonging to particular Williamson sets.
	     */

	    private void calTotalAcidsWillSets() {
		
		Map<String,HashSet<Character>> sets = ConservationSets.williamsonSets();

		Map<String,Integer> setsFreq = new HashMap<String,Integer>();
		
		Set<String> setsKeys = sets.keySet();
		
		Iterator<String> setsKeysItr = setsKeys.iterator();
		
		if (totalFrequency == null) { 
			
			this.calTotalAcidsFrequency();
		}
		
		Map<Character,Integer> totalFreq = totalFrequency;
		
		Set<Character> totalFreqKeys = totalFreq.keySet();
		
		while(setsKeysItr.hasNext()) {
			
			String setsKey = setsKeysItr.next();
		
			Iterator<Character> totalFreqItr = totalFreqKeys.iterator();
		
				while(totalFreqItr.hasNext()) {
					
	                Character totalFreqKey = totalFreqItr.next();
					
					if (sets.get(setsKey).contains(totalFreqKey)) {
						
						Integer count = setsFreq.get(setsKey);
						
						if (count == null) {
							
							setsFreq.put(setsKey, totalFreq.get(totalFreqKey));
						}
						
						else {
							
							setsFreq.put(setsKey, count + totalFreq.get(totalFreqKey));
							
						}
					
					}
				
				}
			
		}
		
		willSetsTotal = setsFreq;
		
		}
	    
	    /**
	     * Calculates a weight of particular sequence according to Vingron-Argos model.
	     * @param seqNr
	     * @return sequence weight
	     */
	    
	    private double weightOfSequenceVingronArgos (int seqNr) {
			
			double weight = 0.0;
			
			for( int i = 0; i < this.numberOfRows(); i++) {
				
				if (i != seqNr) {
					
					weight += ConservationAccessory.percentIdentity(this.getRow(seqNr), this.getRow(i));
					
				}
			
			}
			
			double result  = (1.0 / this.numberOfRows()) * weight ;
			
			return result;
			
		}
		
		/**
		 * Calculates the weight for all the sequences in the alignment.
		 * Weight calculated according to Vingron-Argos model
		 * 
		 */
	    
	    private void weightOfSequencesVingronArgos() {
	    	
	    	vingronArgosWeights = new double[this.numberOfRows()];
	    	
	    	for (int i = 0; i < this.numberOfRows(); i++) {
	    		
	    		vingronArgosWeights[i] = this.weightOfSequenceVingronArgos(i);
	    	}
	    }
	    
	    /**
	     * Gets the values of  Vingron-Argos weights for the whole alignment.
	     * Indices in the weights array correspond to the indices of the sequences in the matrix.  
	     * 
	     * @return an array of weights, indices correspond to sequence numbers
	     */
	    
	    double[] vingronArgosWeights() {
	    	
	    	if (vingronArgosWeights == null) {
	    		
	    		this.weightOfSequencesVingronArgos();
	    		
	    	}
	    	
	    	return vingronArgosWeights;
	    	
	    	}
	    
	    /**
	     * Calculates percent identity for all the sequences in the alignment.
	     * Stores calculated values in a 2D array. The sequence number index the values in the array.
	     * For example percentage identity of sequence 0 and 5 perIden[0][5] 
	     */
	    
	    private void calPercentIdentity() {
	    	
	    	percentIdentity = new double[this.numberOfRows()][this.numberOfRows()];
	    	
	    	int ident = 0;
			
	    		for( int i = 0; i < this.numberOfRows(); i++) {
	    		
	    			for(int j = i + 1; j < this.numberOfRows(); j++){
			
	    				for (int a = 0; a < this.numberOfColumns(); a++) {
				
	    					if (this.getRow(i) == this.getRow(j)) {
					
	    							ident++;
	    					}
			
	    				}
	    				
	    				double result = (double) ident / this.numberOfColumns();
	    				
	    				percentIdentity[i][j] = result;
	    				
	    			}
	    		
	    		}
	    		
	    }
	    
	    /**
	     * Gets percentage identity for the whole alignment.
	     * @return 2D array; the sequences' numbers index the values in the array.
	     */
			
			double[][] getPercentIdentity() {
				
				if(percentIdentity == null) {
					
					this.calPercentIdentity();
					
				}
				
				return percentIdentity;
				
			}
				
			/**
			 * Scary, scary method. Calculates weights according to Voronoi scheme.
			 * @param iter
			 */
		
		void voronoiWeights( int iterNr) {
				
				
				int iterations = iterNr;
				
				Random rgen = new Random();
				
				double[] weights = new double[numberOfRows()];
				
				char[] randSeq = new char[numberOfColumns()];
				
				// sets iterations, don't really know what does it to, but jon set up the iterations to 1000
				
				for ( int i = 0; i < iterations; i++) {
				
					// generates a random sequence, equal in length to the sequences in the alignment
					
					for (int j = 0; j < numberOfColumns(); j++) {
						
						int random = rgen.nextInt(numberOfRows());
						
						randSeq[j] = getMatrixPosition(random, j);
						
					}
					
					// measure the distance between each sequence and a random sequence generated
					// distance measured as percentage identity
					
					double[] distances = new double[numberOfRows()];
					
					double closestValue = 0;
					
					for (int a = 0; a < numberOfRows(); a++) {
						
						distances[a] = 1.0 - ConservationAccessory.percentIdentity(getRow(a), randSeq);
						
						if(distances[a] < closestValue) {
							
							closestValue = distances[a];
							
						}
						
					}
					
					// collect all the sequences with the closest distance 
					
					List<Integer> closestSeqs = new ArrayList<Integer>();
					
					for (int b = 0; b < distances.length; b++) {
						
						double dis = distances[b];
						
						if ( dis == closestValue) {
							
							closestSeqs.add(b);
						}
					}
					
					// increase by one the weight of the closest sequence
					
					double increase = 1.0 / closestSeqs.size();
					
					for (int c = 0; c < closestSeqs.size(); c++ ) {
						
						int cs = closestSeqs.get(c);
						
						weights[cs] += increase;
						
					}
					
					// repeat iterations times
					
				}
				
				// normalize weights so they sum up to N
				
				double weightSum = 0.0;
				
				for (int d = 0; d < weights.length; d++) {
					
					weightSum += weights[d];
					
				}
				
				double scaleFactor = weightSum / numberOfRows();
				
				for (int e = 0; e < weights.length; e++) {
					
					weights[e] = weights[e] + scaleFactor;
					
				}
				
				voronoiWeighths = weights;
				
			}
			
			/** 
			 * Gets Voronoi weights.
			 * @param iterNr
			 * @return array containing voronoi weights, indices correspond to sequence indices in matrix
			 */
			double[] getVoronoiWeights(int iterNr) {
				
				if (voronoiWeighths == null) {
					
					this.voronoiWeights(iterNr);
				}
			
			return voronoiWeighths;
			
			}
		
	
		
		
	}
			
			
			
			
			
			





