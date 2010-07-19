package compbio.examples;

import java.util.*;
public class SmallestElement {
	
	public static void main (String[] args) {
		
		Integer one = new Integer(1);
		Integer two = new Integer(2);
		Integer three = new Integer(3);
		
		Map<String,Integer> map = new HashMap<String, Integer>();
		
		map.put("one", one);
		map.put("two",two );
		map.put("three", three);
		
		Collection<Integer> c = map.values();
		
		int min = Collections.min(c);
		
		System.out.println(min);
	}

}
