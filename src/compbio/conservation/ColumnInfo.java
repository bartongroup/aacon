package compbio.conservation;

import java.io.*;

public class ColumnInfo {
	
	private final String group;
	
	private final double consScore;
	
	private final String properties;
	
	public ColumnInfo (String group, double score, String properties) {
		
		this.group = group;
		
		this.consScore = score;
		
		this.properties = properties;
		
	}
	
	String getGroup() {
		
		return group;
	}
	
	double getScore() {
		
		return consScore;
		
	}
	
	String getProperties() {
		
		return properties;
	}
	
	String getStatus(Method method) {
		
		String stat = ConservationStatus.stringReps(ConservationStatus.getStatus(consScore, method));
			
			return stat;
	}
	
	void printInfo(Method method, PrintWriter print) {
		
		int width1 = 20;
		
		int width2 = 20;
		
		String format1 = "%-" + width1 + "s";
		
		String format2 = "%-" + width2 + "s";
		
		print.printf(format1, group);
		
		print.printf(format2, getStatus(method));
		
		print.println(properties);
	}

}
