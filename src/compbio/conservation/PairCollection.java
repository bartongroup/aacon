package compbio.conservation;

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

		}




