package compbio.conservation;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;

public class ShannonEnthropy {

	
	static double ShannonLog2(final Map<Character,Integer> map, double nrSequences) {
		
		if (map == null) {
			throw new IllegalArgumentException("Map must not be null");
			}
		
		assert !map.isEmpty() : "Shannon has been fed an empty map";
		assert nrSequences > 0;
		
        double sum = 0;
    	
    	Collection<Character> keys = map.keySet();
    	
    	Iterator<Character> itr = keys.iterator();
    	
    	while(itr.hasNext()) {
    		
    		if (itr.next() != '-') {
    			
    			double value = map.get(itr.next()) / nrSequences;
    			
    			sum = sum + (value * (Math.log(value)/Math.log(2.0)));
    			
    		}
    			
    	}
    		
    	
    	assert sum > 0 : "Shannon has been fed an empty column";
    	
    	return  -sum;
    	
    	
	}

static double ShannonLn(final Map<Character,Integer> map, int nrSequences) {
		
		if (map == null) {
			throw new IllegalArgumentException("Map must not be null");
			}
		
		assert !map.isEmpty() : "Shannon has been fed an empty map";
		assert nrSequences > 0; 
		
        double sum = 0;
    	
    	Collection<Character> keys = map.keySet();
    	
    	Iterator<Character> itr = keys.iterator();
    	
    	while(itr.hasNext()) {
    		
    		if (itr.next() != '-') {
    			
    			double value = map.get(itr.next()) / (double) nrSequences;
    			
    			sum = sum + (value * Math.log(value));
    			
    		}
    			
    	}
    		
    	
    	assert sum > 0 : "Shannon enthropy has been fed an empty column";
    	
    	return  -sum;
    	
    	
	}

	
	
}
