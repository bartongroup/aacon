package compbio.conservation;

import java.util.*;
import java.io.*;
import compbio.util.*;

public class _WrapperClass2 {
	
	
	
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
		
		
	AminoAcidMatrix matrix = new AminoAcidMatrix(fastaSeqs);
		
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
	
	String fileName = "results.txt";
		
	
	int nrOfCols = matrix.numberOfColumns();
	
	int winWidth = 7;
	
	boolean kabat = false; 
	
	boolean jores = false;  
	
	boolean schneider = false; 
	
	boolean shenkin = false; 
	
	boolean gerstein = false;
	
	boolean taylorNoGaps = false;
	
	boolean taylorGaps = false;
	
	boolean zvelibil = false;
	
	boolean karlin = false;
	
	boolean armon = false;
	
	boolean thompson = false;
	
	final boolean lancet  = false;
	
	boolean mirny = true;
	
	boolean williamson = false;
	
	boolean landgraf = false;
	
	boolean sander = false;
	
	boolean valdar = false;
	
	List<String> scores = new ArrayList<String>();
	
	if(kabat) {
		
		scores.add("Kabat");
	
	}
		
	if(jores) {
		
		scores.add("Jores");
		
	}
	
	if(schneider) {
		
		scores.add("Schneider");
		
	}
	
	if(shenkin) {
		
		scores.add("Shenkin");
		
	}
	
	if(gerstein) {
		
		scores.add("Gerstein");
		
	}
	
	if(taylorNoGaps) {
		
		scores.add("TaylorNoGaps");
		
	}
	if(taylorGaps) {
		
		scores.add("TaylorGaps");
		
	}
	
	if(zvelibil) {
		
		scores.add("Zvelibil");
		
	}
		
	if(karlin) {
		
		scores.add("Karlin");
		
	}
	
	if(armon) {
		
		scores.add("Armon");
		
	}
	
	if(thompson) {
		
		scores.add("Tompson");
		
	}
	
	if(lancet) {
		
		scores.add("Lancet");
		
	}
		
	if(mirny) {
		
		scores.add("Mirny");
		
	}
	
	if(williamson) {
		
		scores.add("Williamson");
		
	}
		
	if(landgraf) {
		
		scores.add("Landgraf");
	}
	
	if (sander) {
		
		scores.add("Sander");
		
	}
	
	if(valdar) {
		
		scores.add("Valdar");
	
	}
	
	PrintWriter print = null;
	
	try {
		
		print = new PrintWriter( new BufferedWriter (new FileWriter(fileName)));
	}
	
	catch(IOException ex) {
		
		System.out.println("Problem writing" + fileName);
		
	}
	
	int colLen = matrix.numberOfRows();
	
	int targetWidth = colLen + 10;
	
	String formatString = "%-"+targetWidth+"s";
	
	print.printf(formatString, "Columns");
	
	for( int l = 0; l < scores.size(); l++){
		
		print.printf("%-30s", scores.get(l));
	}
	
	print.println("");
	
	_ColumnCollectionShortList cols = new _ColumnCollectionShortList(matrix,winWidth);
	
	//ColumnCollectionList colsRem = null;
	
	cols.calculationInitializtion(kabat, jores, schneider, shenkin, gerstein, taylorNoGaps, taylorGaps, zvelibil, karlin, armon, thompson, lancet, mirny, williamson, landgraf, sander, valdar);
	
	for(int z = 0; z < winWidth; z++) {
		
		print.printf(formatString, z+ " " + new String(cols.getColumnCollection().get(z).getColumnAcids()));
		
		if(kabat) {
			
			print.printf("%-30.6f", cols.getKabat().get(z));
	
		}
		
		if(jores) {
		
			print.printf("%-30.6f", cols.getJores().get(z));
		
		}
	
		if(schneider) {
		
			print.printf("%-30.6f", cols.getSchneider().get(z));
		
		}
	
		if(shenkin) {
		
			print.printf("%-30.6f", cols.getShenkin().get(z));
		}
	
		if(gerstein) {
		
			print.printf("%-30.6f", cols.getGerstein().get(z));
		}
	
		if(taylorNoGaps) {
		
			print.printf("%-30d", cols.getTaylorNoGaps().get(z));
		}
		
		if(taylorGaps) {
		
			print.printf("%-30d", cols.getTaylorGaps().get(z));
		
		}
	
		if(zvelibil) {
		
			print.printf("%-30d", cols.getZvelibil().get(z));
		
		}
		
		if(karlin) {
		
			print.printf("%-30.6f", cols.getKarlin().get(z));
		
		}
	
		if(armon) {
		
			print.printf("%-30.6f", cols.getArmon().get(z));
		}
	
		if(thompson) {
		
			print.printf("%-30.6f", cols.getThompson().get(z));
		}
	
		if(lancet) {
		
			print.printf("%-30.6f", cols.getLancet().get(z));
	
		}
		
		if(mirny) {
		
			print.printf("%-30.6f", cols.getMirny().get(z));
		
		}
	
		if(williamson) {
		
			print.printf("%-30.6f", cols.getWilliamson().get(z));
	
		}
		
		if(landgraf) {
		
			print.printf("%-30.6f", cols.getLandgraf().get(z));
		
		}
	
		if (sander) {
		
			print.printf("%-30.6f", cols.getSander().get(z));
		
		}
	
		if(valdar) {
		
			print.printf("%-30.6f%", cols.getKabat().get(z));
		}
		
		print.println("");
	
	}
	
	//int rem = nrOfCols%winWidth;
	
	for (int k = winWidth; k < nrOfCols; k++ ) {
		
		cols.changeColumns(k);
		
		cols.calculationUpdate(kabat, jores, schneider, shenkin, gerstein, taylorNoGaps, taylorGaps, zvelibil, karlin, armon, thompson, lancet, mirny, williamson, landgraf, sander, valdar);
		
		print.printf(formatString, k + " "+ new String(cols.getNewestColumn().getColumnAcids()));
		
			if(kabat) {
			
				print.printf("%-30.6f", cols.getKabat().get(winWidth - 1));
		
			}
			
			if(jores) {
			
				print.printf("%-30.6f", cols.getJores().get(winWidth - 1));
			
			}
		
			if(schneider) {
			
				print.printf("%-30.6f", cols.getSchneider().get(winWidth - 1));
			
			}
		
			if(shenkin) {
			
				print.printf("%-30.6f", cols.getShenkin().get(winWidth - 1));
			}
		
			if(gerstein) {
			
				print.printf("%-30.6f", cols.getGerstein().get(winWidth - 1));
			}
		
			if(taylorNoGaps) {
			
				print.printf("%-30d", cols.getTaylorNoGaps().get(winWidth - 1));
			}
			
			if(taylorGaps) {
			
				print.printf("%-30d", cols.getTaylorGaps().get(winWidth - 1));
			
			}
		
			if(zvelibil) {
			
				print.printf("%-30d", cols.getZvelibil().get(winWidth - 1));
			
			}
			
			if(karlin) {
			
				print.printf("%-30.6f", cols.getKarlin().get(winWidth - 1));
			
			}
		
			if(armon) {
			
				print.printf("%-30.6f", cols.getArmon().get(winWidth - 1));
			}
		
			if(thompson) {
			
				print.printf("%-30.6f", cols.getThompson().get(winWidth - 1));
			}
		
			if(lancet) {
			
				print.printf("%-30.6f", cols.getLancet().get(winWidth - 1));
		
			}
			
			if(mirny) {
			
				print.printf("%-30.6f", cols.getMirny().get(winWidth - 1));
			
			}
		
			if(williamson) {
			
				print.printf("%-30.6f", cols.getWilliamson().get(winWidth - 1));
		
			}
			
			if(landgraf) {
			
				print.printf("%-30.6f", cols.getLandgraf().get(winWidth - 1));
			
			}
		
			if (sander) {
			
				print.printf("%-30.6f", cols.getSander().get(winWidth - 1));
			
			}
		
			if(valdar) {
			
				print.printf("%-30.6f", cols.getKabat().get(winWidth - 1));
		}
			
			print.println("");
		
	}
	
	print.close();
	
	}
	
}
	
	