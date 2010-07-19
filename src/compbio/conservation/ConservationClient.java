package compbio.conservation;

import java.io.*;
import java.util.*;

import compbio.util.FastaSequence;
import compbio.util.SequenceUtil;

class ConservationClient {
	
	//static enum Method { kabatScore, joresScore, schneiderScore, shenkinScore, gersteinScore, SmallestTaylorSetGaps, SmallestTaylorSetNoGaps, zvelibilScore, karlinScore, armonScore, thompsonScore, notLancetScore, mirnyScore, williamsonScore, landgrafScore, sanderScore, valdarScore;
	
		//static Method getMethod(String meth) {
		
			//meth = meth.trim().toLowerCase();
		
			//if(meth.equalsIgnoreCase(kabatScore.toString())) {
			
				//return kabatScore;
				
			//}
			
			//if(meth.equalsIgnoreCase(joresScore.toString())) {
				
				//return joresScore;
				
			///}
			
			//if(meth.equalsIgnoreCase(schneiderScore.toString())) {
				
				//return schneiderScore;
				
			//}
			
			//if(meth.equalsIgnoreCase(shenkinScore.toString())) {
				
				//return shenkinScore;
				
			//}
			
			//if(meth.equalsIgnoreCase(gersteinScore.toString())) {
				
				//return gersteinScore;
				
			//}
			
			//if(meth.equalsIgnoreCase(SmallestTaylorSetGaps.toString())) {
				
				//return SmallestTaylorSetGaps;
				
			//}
			
			//if(meth.equalsIgnoreCase(SmallestTaylorSetNoGaps.toString())) {
				
				//return SmallestTaylorSetNoGaps;
				
			//}
			
			//if(meth.equalsIgnoreCase(zvelibilScore.toString())) {
				
				//return zvelibilScore;
				
			//}
			
			//if(meth.equalsIgnoreCase(karlinScore.toString())) {
				
				//return karlinScore;
				
			//}
			
			//if(meth.equalsIgnoreCase(armonScore.toString())) {
				
				//return armonScore;
				
			//}
			
			//if(meth.equalsIgnoreCase(thompsonScore.toString())) {
				
				//return thompsonScore;
				
			//}
			
			//if(meth.equalsIgnoreCase(notLancetScore.toString())) {
				
				//return notLancetScore;
				
			//}
			
			//if(meth.equalsIgnoreCase(mirnyScore.toString())) {
				
				//return mirnyScore;
				
			//}
			
			//if(meth.equalsIgnoreCase(williamsonScore.toString())) {
				
				//return williamsonScore;
				
			//}
			
			//if(meth.equalsIgnoreCase(landgrafScore.toString())) {
				
				//return landgrafScore;
				
			//}
			
			//if(meth.equalsIgnoreCase(sanderScore.toString())) {
				
				//return sanderScore;
				
			//}
			
			//if(meth.equalsIgnoreCase(valdarScore.toString())) {
				
				//return valdarScore;
				
			//}
			
			//return null;
	
		//}
	
	//}
	
	private final Map<Method, double[]> scores = new EnumMap<Method, double[]>(Method.class);
	
	final static String pseparator = "=";
	
	final static String methodKey = "-m";
	
	final static String normalizationKey = "-n";
	
	final static String formatKey = "-f";
	
	final static String inputKey = "-i";
	
	final static String outputKey = "-o";
	
	static String[] getMethodNames(String[] cmd) {
		
		for (int i = 0; i < cmd.length; i++) {
			
			String meths = cmd[i];
			
			if(meths.trim().toLowerCase().startsWith(methodKey + pseparator)) {
				
				return meths.substring(meths.indexOf(pseparator) + 1).split(",");
			}
		}
		
		return null;
		
	}
	
	static boolean getNormalize(String[] cmd) {
		
		for (int i = 0; i < cmd.length; i++) {
			
			String norm = cmd[i];
			
			if(norm.trim().toLowerCase().equals(normalizationKey + pseparator)) {
				
				return true;
			}
		}
		
		return false;
		
	}
	
	static String getFormat(String[] cmd) {
		
		for (int i = 0; i < cmd.length; i++) {
			
			String form = cmd[i];
			
			if(form.trim().toLowerCase().startsWith(formatKey + pseparator)) {
				
				return form.substring(form.indexOf(pseparator) + 1);
			}
		}
		
		return null;
	}
	
	static String getOutputFilePath(String[] cmd) {
		
		for (int i = 0; i < cmd.length; i++) {
			
			String name = cmd[i];
			
			if(name.trim().toLowerCase().startsWith(outputKey + pseparator)) {
				
				return name.substring(name.indexOf(pseparator) + 1);
			}
		}
		
		return null;
	}
	
	static String getInputFilePath(String[] cmd) {
		
		for (int i = 0; i < cmd.length; i++) {
			
			String name = cmd[i];
			
			if(name.trim().toLowerCase().startsWith(inputKey + pseparator)) {
				
				return name.substring(name.indexOf(pseparator) + 1);
			}
		}
		
		return null;
	}
	
	
	
	// don't know how to set up a path, so far path read from the command line;
	// but probably it is not how it is normally done
	// waiting for suggestions
	// formats not yet decided
	
	double[] calculateMethod(Method method, AminoAcidMatrix matrix, boolean normalize) {
		
		double[] result = null;
		
		if (method == null) {
			
			System.out.println("Method not suppoted");
			
			return result;
		}
	
		else {
			
			ConservationScores scores = new ConservationScores(matrix);
			
			if(method.equals(Method.kabatScore)) {
				
				result =  scores.kabatScore(normalize);
			}
			
			if(method.equals(Method.joresScore)) {
				
				result =  scores.joresScore(normalize);
			}
			
			if(method.equals(Method.schneiderScore)) {
				
				result =  scores.schneiderScore(normalize);
			}
			
			if(method.equals(Method.shenkinScore)) {
				
				result =  scores.shenkinScore(normalize);
			}
			
			if(method.equals(Method.gersteinScore)) {
				
				result =  scores.gersteinScore(normalize);
			}
			
			if(method.equals(Method.SmallestTaylorSetGaps)) {
				
				result =  scores.SmallestTaylorSetGaps(normalize);
			}
			
			if(method.equals(Method.SmallestTaylorSetNoGaps)) {
				
				result=  scores.SmallestTaylorSetNoGaps(normalize);
			}
			
			if(method.equals(Method.zvelibilScore)) {
				
				result =  scores.zvelibilScore(normalize);
			}
			
			if(method.equals(Method.karlinScore)) {
				
				result =  scores.karlinScore(normalize);
			}
			
			if(method.equals(Method.armonScore)) {
				
				result =  scores.armonScore(normalize);
			}
			
			if(method.equals(Method.thompsonScore)) {
				
				result =  scores.thompsonScore(normalize);
			}
			
			if(method.equals(Method.notLancetScore)) {
				
				result =  scores.notLancetScore(normalize);
			}
			
			if(method.equals(Method.mirnyScore)) {
				
				result =  scores.mirnyScore(normalize);
			}
			
			if(method.equals(Method.williamsonScore)) {
				
				result =  scores.williamsonScore(normalize);
			}
			
			if(method.equals(Method.landgrafScore)) {
				
				result =  scores.landgrafScore(normalize);
			}
			
			if(method.equals(Method.sanderScore)) {
				
				result =  scores.sanderScore(normalize);
			}
			
			if(method.equals(Method.valdarScore)) {
				
				result =  scores.valdarScore(normalize);
			}
			
			return result;
			
		}
	
	}
	
	ConservationClient(String[] cmd) {
		
		String[] methods = getMethodNames(cmd);
		
		if(methods == null) {
			
			System.out.println("Methods not provided");
			
		}
		
		String format = getFormat(cmd);
		
		if (format == null) {
			
			System.out.println("Format not provided");
			
		}
		
		String inFilePath = getInputFilePath(cmd);
		
		if (inFilePath == null) {
			
			System.out.println("Input file path not provided.");
			
		}
		
		String outFilePath = getOutputFilePath(cmd);
		
		if (outFilePath == null) {
			
			System.out.println("Output file path not provided.");
			
		}
		
		boolean normalize = getNormalize(cmd);
		
		InputStream inStr = null;
		
		List<FastaSequence> fastaSeqs = null;
		
		try {
			
			inStr = new FileInputStream(inFilePath);
			
		}
		
		catch (IOException e) {
			
			System.out.println("Can not find file");
		
		}
		
		try {
			
			fastaSeqs = SequenceUtil.readFasta(inStr);
		}
		
		catch (IOException e) {
			
			System.out.println("Sth wrong with reading the file");
		}
			
		
		AminoAcidMatrix alignment = new AminoAcidMatrix(fastaSeqs);
		
		for (int i = 0; i < methods.length; i++) {
			
			Method meth = Method.getMethod(methods[i]);
			
			scores.put(meth, this.calculateMethod(meth, alignment, normalize));
			
		}
		
		ConservationFormatter.formatResults(scores);
		
		}
	
	public static void main(String[] args) {
		
		if(args == null) {
			
			System.out.println ("No parameters were suppled");
		}
		
		if(args.length < 6) {
			
			System.out.println("Method names, normalizetion status, output format, input and output file paths are required.");
		}
		
		ConservationClient cons = new ConservationClient(args);
	}

}
