package compbio.conservation;

import java.util.*;
import java.io.*;
import compbio.util.*;

public class _WrapperClass {
	
	static void printListToFileDoub(List<Double> listName, String fileName, String scoreName) {
		
		PrintWriter print = null;
		
		List<Double> myList = listName;
		
		try {
			
			print = new PrintWriter( new BufferedWriter (new FileWriter(fileName)));
		}
		
		catch(IOException e) {
			
			System.out.println("Problem writing" + fileName);
			
		}
		
		print.println(scoreName);
		
		for (int i = 0; i < myList.size(); i++ ) {
			
			print.println("Column " + i + ": " + myList.get(i));
			
			}
		
		print.close();
	}
	
	static void printListToFileInt(List<Integer> listName, String fileName, String scoreName) {
		
		PrintWriter print = null;
		
		List<Integer> myList = listName;
		
		try {
			
			print = new PrintWriter( new BufferedWriter (new FileWriter(fileName)));
		}
		
		catch(IOException e) {
			
			System.out.println("Problem writing" + fileName);
			
		}
		
		print.println(scoreName);
		
		for (int i = 0; i < myList.size(); i++ ) {
			
			print.println("Column " + i + ": " + myList.get(i));
			
			}
		
		print.close();
		
	}
	
	
	public static void main (String[] args) {

	
	String filePath  = "/homes/agolicz/alignments/alignment1";
	
	InputStream inStr = null;
	
	List<FastaSequence> fastaSeqs = null;
	
	try {
		
		inStr = new FileInputStream(filePath);
		
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
		
		
	AminoAcidMatrix matrix = new AminoAcidMatrix(fastaSeqs, null);
		
		//char a = 'D';
		//char b = 'D';
		//char c = 'D';
		//char d = 'D';
		//char e = 'D';
		//char f = 'E';
		//char g = 'E';
		//char h = 'E';
		//char i = 'F';
		//char j = '-';
	
	//AminoAcidMatrix matrix = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i, j);
		
	
	int nrOfCols = matrix.numberOfColumns();
	
	int winWidth = 7;
	
	boolean kabat = true; 
	
	boolean jores = true;  
	
	boolean schneider = true; 
	
	boolean shenkin = true; 
	
	boolean gerstein = true;
	
	boolean taylorNoGaps = true;
	
	boolean taylorGaps = true;
	
	boolean zvelibil = true;
	
	boolean karlin = true;
	
	boolean armon = true;
	
	boolean thompson = true;
	
	boolean lancet = true;
	
	boolean mirny = true;
	
	boolean williamson = true;
	
	boolean landgraf = true;
	
	boolean sander = true;
	
	boolean valdar = true;
	
	
	_ColumnCollectionList cols = new _ColumnCollectionList(matrix,winWidth);
	
	_ColumnCollectionList colsRem = null;
	
	cols.calculationInitializtion(kabat, jores, schneider, shenkin, gerstein, taylorNoGaps, taylorGaps, zvelibil, karlin, armon, thompson, lancet, mirny, williamson, landgraf, sander, valdar);
	
	int rem = nrOfCols%winWidth;
	
	for (int k = winWidth; k < nrOfCols - rem; k = k + winWidth ) {
		
		cols.changeColumns(k);
		
		cols.calculationUpdate(kabat, jores, schneider, shenkin, gerstein, taylorNoGaps, taylorGaps, zvelibil, karlin, armon, thompson, lancet, mirny, williamson, landgraf, sander, valdar);
		
	}
	
	if (rem == 0) {

		if(kabat) {
			
			List<Double> result = cols.getKabat();
			
			assert !result.isEmpty();
			
			printListToFileDoub(result,"kabat.txt", "Kabat" );
		
		}
			
		if(jores) {
			
			List<Double> result = cols.getJores();
			
			printListToFileDoub(result,"jores.txt", "Jores" );
			
		}
		
		if(schneider) {
			
			List<Double> result = cols.getSchneider();
			
			printListToFileDoub(result,"schneider.txt", "Schneider" );
			
		}
		
		if(shenkin) {
			
			List<Double> result = cols.getShenkin();
			
			printListToFileDoub(result,"shenkin.txt", "Shenkin" );
			
		}
		
		if(gerstein) {
			
			List<Double> result = cols.getGerstein();
			
			printListToFileDoub(result,"gerstein.txt", "Gerstein" );
			
		}
		
		if(taylorNoGaps) {
			
			List<Integer> result = cols.getTaylorNoGaps();
			
			printListToFileInt(result,"tylorNoGaps.txt", "tylorNoGaps" );
			
		}
		if(taylorGaps) {
			
			List<Integer> result = cols.getTaylorGaps();
			
			printListToFileInt(result,"talorGaps.txt", "tylorGaps" );
			
		}
		
		if(zvelibil) {
			
			List<Integer> result = cols.getZvelibil();
			
			printListToFileInt(result,"zvelibil.txt", "Zvelibil" );
			
		}
			
		if(karlin) {
			
			List<Double> result = cols.getKarlin();
			
			printListToFileDoub(result,"karlin.txt", "Kalin" );
			
		}
		
		if(armon) {
			
			List<Double> result = cols.getArmon();
			
			printListToFileDoub(result,"armon.txt", "Armon" );
			
		}
		
		if(thompson) {
			
			List<Double> result = cols.getThompson();
			
			printListToFileDoub(result,"thompson.txt", "Thompson" );
			
		}
		
		if(lancet) {
			
			List<Double> result = cols.getLancet();
			
			printListToFileDoub(result,"lancet.txt", "Lancet" );
			
		}
			
		if(mirny) {
			
			List<Double> result = cols.getMirny();
			
			printListToFileDoub(result,"mirny.txt", "Mirny" );
			
		}
		
		if(williamson) {
			
			List<Double> result = cols.getWilliamson();
			
			printListToFileDoub(result,"williamson.txt", "Willimason" );
			
		}
			
		if(landgraf) {
			
			List<Double> result = cols.getLandgraf();
			
			printListToFileDoub(result,"landgraf.txt", "Landgraf" );
		}
		
		if (sander) {
			
			List<Double> result = cols.getSander();
			
			printListToFileDoub(result,"sander.txt", "Sander" );
			
		}
		
		if(valdar) {
			
			List<Double> result = cols.getValdar();
			
			printListToFileDoub(result,"valdar.txt", "Valdar" );
		
		}
		
	}
	
	else {
		
		colsRem = new _ColumnCollectionList(matrix, rem);
		
		colsRem.calculationInitializtion(kabat, jores, schneider, shenkin, gerstein, taylorNoGaps, taylorGaps, zvelibil, karlin, armon, thompson, lancet, mirny, williamson, landgraf, sander, valdar);

		if(kabat) {
			
			Collection<Double> remResult = colsRem.getKabat();
			
			List<Double> result = cols.getKabat();
			
			result.addAll(remResult);
			
			printListToFileDoub(result,"kabat.txt", "Kabat" );
		
		}
			
		if(jores) {
			
			Collection<Double> remResult = colsRem.getJores();
			
			List<Double> result = cols.getJores();
			
			result.addAll(remResult);
			
			printListToFileDoub(result,"jores.txt", "Jores" );
			
		}
		
		if(schneider) {
			
			Collection<Double> remResult = colsRem.getSchneider();
			
			List<Double> result = cols.getSchneider();
			
			result.addAll(remResult);
			
			printListToFileDoub(result,"schneider.txt", "Schneider" );
			
		}
		
		if(shenkin) {
			
			Collection<Double> remResult = colsRem.getShenkin();
			
			List<Double> result = cols.getShenkin();
			
			result.addAll(remResult);
			
			printListToFileDoub(result,"shenkin.txt", "Shenkin" );
			
		}
		
		if(gerstein) {
			
			Collection<Double> remResult = colsRem.getGerstein();
			
			List<Double> result = cols.getGerstein();
			
			result.addAll(remResult);
			
			printListToFileDoub(result,"gerstein.txt", "Gerstein" );
			
		}
		
		if(taylorNoGaps) {
			
			Collection<Integer> remResult = colsRem.getTaylorNoGaps();
			
			List<Integer> result = cols.getTaylorNoGaps();
			
			result.addAll(remResult);
			
			printListToFileInt(result,"tylorNoGaps.txt", "tylorNoGaps" );
			
		}
		if(taylorGaps) {
			
			Collection<Integer> remResult = colsRem.getTaylorGaps();
			
			List<Integer> result = cols.getTaylorGaps();
			
			result.addAll(remResult);
			
			printListToFileInt(result,"talorGaps.txt", "tylorGaps" );
			
		}
		
		if(zvelibil) {
			
			Collection<Integer> remResult = colsRem.getZvelibil();
			
			List<Integer> result = cols.getZvelibil();
			
			result.addAll(remResult);
			
			printListToFileInt(result,"kabat.txt", "Kabat" );
			
		}
			
		if(karlin) {
			
			Collection<Double> remResult = colsRem.getKarlin();
			
			List<Double> result = cols.getKarlin();
			
			result.addAll(remResult);
			
			printListToFileDoub(result,"karlin.txt", "Kalin" );
			
		}
		
		if(armon) {
			
			Collection<Double> remResult = colsRem.getArmon();
			
			List<Double> result = cols.getArmon();
			
			result.addAll(remResult);
			
			printListToFileDoub(result,"armon.txt", "Armon" );
			
		}
		
		if(thompson) {
			
			Collection<Double> remResult = colsRem.getThompson();
			
			List<Double> result = cols.getThompson();
			
			result.addAll(remResult);
			
			printListToFileDoub(result,"thompson.txt", "Thompson" );
			
		}
		
		if(lancet) {
			
			Collection<Double> remResult = colsRem.getLancet();
			
			List<Double> result = cols.getLancet();
			
			result.addAll(remResult);
			
			printListToFileDoub(result,"lancet.txt", "Lancet" );
			
		}
			
		if(mirny) {
			
			Collection<Double> remResult = colsRem.getMirny();
			
			List<Double> result = cols.getMirny();
			
			result.addAll(remResult);
			
			printListToFileDoub(result,"mirny.txt", "Mirny" );
			
		}
		
		if(williamson) {
			
			Collection<Double> remResult = colsRem.getWilliamson();
			
			List<Double> result = cols.getWilliamson();
			
			result.addAll(remResult);
			
			printListToFileDoub(result,"williamson.txt", "Willimason" );
			
		}
			
		if(landgraf) {
			
			Collection<Double> remResult = colsRem.getLandgraf();
			
			List<Double> result = cols.getLandgraf();
			
			result.addAll(remResult);
			
			printListToFileDoub(result,"landgraf.txt", "Landgraf" );
		}
		
		if (sander) {
			
			Collection<Double> remResult = colsRem.getSander();
			
			List<Double> result = cols.getSander();
			
			result.addAll(remResult);
			
			printListToFileDoub(result,"sander.txt", "Sander" );
			
		}
		
		if(valdar) {
			
			Collection<Double> remResult = colsRem.getValdar();
			
			List<Double> result = cols.getValdar();
			
			result.addAll(remResult);
			
			printListToFileDoub(result,"valdar.txt", "Valdar" );
		
		}
		
	}
	
	}

}
	
	
	
	