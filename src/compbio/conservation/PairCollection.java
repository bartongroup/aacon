package compbio.conservation;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;


// some of the methods from column collection got dumped here in case they will be needed, not likely though

	public class PairCollection {
		
		private int nrPairs;

		private int mostFreqPair;

		public PairCollection(int p, int m) {
			
			nrPairs = p;
			
			mostFreqPair = m;
		}

		public int getNrPairs() {
			
			return nrPairs;
			
		}

		public int getMostFrequent() {
			
			return mostFreqPair;
			
		}

		

	Map<Character,Integer> TotalAcids() {
		
		Map<Character,Integer>  freq = new HashMap<Character,Integer>();
		
		
		return freq;
		
	}

	// calculate the talat nuber of aa in the alignament that belon to a particular set

	Map<String,Integer> TotalAcidsWillSets() {
		
		Map<String,HashSet<Character>> sets = ConservationAccessory.williamsonSets();

		Map<String,Integer> setsFreq = new HashMap<String,Integer>();
		
		Set<String> setsKeys = sets.keySet();
		
		Iterator<String> setsKeysItr = setsKeys.iterator();
		
		Map<Character,Integer> totalFreq = this.TotalAcids();
		
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
		
		return setsFreq;
		
		}
	}




