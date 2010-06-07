package compbio.conservation;

//creates a matrix of aa in multiple alignment
//gets fasta sequences and puts them into matrix 
//might have to check if all sequences are equal length(ask)

import java.util.List;
import compbio.util.SequenceUtil;
import compbio.util.FastaSequence;
import java.io.InputStream;
import java.io.IOException;

public class AminoAcidMatrix {
	
	private char[][] matrix;
	
// this constructor will be use for testign solely
// it lest the user enter aa into matrix manually 
	
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
	}
	
		
	public AminoAcidMatrix(InputStream inStream){
		
	   List<FastaSequence> seqs = null;
	   
	   try {
	          seqs = SequenceUtil.readFasta(inStream);
	       }
	   catch (IOException e) 
	       {
	          System.out.println("Can not read input Stream");
	       }
	          
	          int sequenceNr = seqs.size();
	           
	          FastaSequence seq = seqs.get(0);
	     
	          String firstSequence = seq.getSequence();

	          int sequenceLength = firstSequence.length();
	       
	          matrix = new char[sequenceNr][sequenceLength];

	          for( int i = 0; i < sequenceNr; i++) {

	                FastaSequence s = seqs.get(i);
	         
	                char[] sequenceChars = s.getSequence().toCharArray();
	                   
			              for ( int j = 0; j < sequenceLength; i++) {
				
	                                 matrix[i][j] = sequenceChars[j];
	                        
	                      }
	          }
	}           

	public int numberOfColumns() {

	int nrColumns = matrix[0].length;

	return nrColumns;

	}

	public int numberOfRows() {

	int nrRows = matrix.length;

	return nrRows;

	}

	public char[][] getMatrix() {

	return matrix;

	}

	public char getMatrixPosition(int row, int column) {

	char position = matrix[row][column];

	return position;

	}
	       
}
