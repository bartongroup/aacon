package compbio.conservation;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.log4j.Logger;

//this is class Column, the objects created represent columns of an alignment
//column collection will be created
//columns are than used in various calculations
//the idea is that all the columns objects will be stored in an array
//that array will be created as instantiation of ColumnCollection class
//will see what's gonna happen

public class AminoAcidColumn {

    private final static Logger log = Logger.getLogger(AminoAcidColumn.class);

    private char[] columnArr;
    private AminoAcidOccurance[] aaOcc;

    // TO BE DELETED, THIS CONSTRUCTOR USED FOR TESTS ONLY

    public AminoAcidColumn(char a, char b, char c, char d, char e, char f,
	    char g, char h, char i, char j) {

	columnArr = new char[10];

	columnArr[0] = a;
	columnArr[1] = b;
	columnArr[2] = c;
	columnArr[3] = d;
	columnArr[4] = e;
	columnArr[5] = f;
	columnArr[6] = g;
	columnArr[7] = h;
	columnArr[8] = i;
	columnArr[9] = j;

	createOccurance();
	countOccurance();

    }

    // constructor, 

    public AminoAcidColumn(AminoAcidMatrix m, int c) {

	int colLength = m.numberOfRows();

	StringBuilder s = new StringBuilder();

	columnArr = new char[colLength];

	// creates a string 

	for (int i = 0; i < colLength; i++) {

	    char aacid = m.getMatrixPosition(i, c);
	    s = s.append(aacid);
	}

	String columnStr = s.toString().toUpperCase();

	columnArr = columnStr.toCharArray();

	// creates occurance objects

	this.createOccurance();

	// counts occurance

	this.countOccurance();

    }

    //returns column length

    public int length() {

	int len = columnArr.length;

	return len;

    }

    //returns occurance

    public AminoAcidOccurance[] getOccuranceArray() {

	return aaOcc;

    }

    // returns Occurance (as integer) of a particular aa we are looking for 

    public int getOccurance(char a) {

	int o = 0;

	for (int i = 0; i < aaOcc.length; i++) {

	    if (aaOcc[i].getId() == a)
		o = aaOcc[i].getOccurance();

	}

	return o;
    }

    // checks whether column is empty(consists of gaps only)

    // FIXME check against fre
    public boolean isEmpty() {

	int nrOfGaps = 0;

	for (int i = 0; i < aaOcc.length; i++) {

	    if (aaOcc[i].getId() == '-')
		nrOfGaps = aaOcc[i].getOccurance();
	}

	if (nrOfGaps == columnArr.length)
	    return true;

	else
	    return false;

    }

    // checks whether there are gaps present in the column,
    //returns true if there are, false if there aren't any

    public boolean gapsInColumn() {

	int nrOfGaps = 0;

	for (int i = 0; i < aaOcc.length; i++) {

	    if (aaOcc[i].getId() == '-')
		nrOfGaps = aaOcc[i].getOccurance();
	}

	if (nrOfGaps > 0)
	    return true;

	else
	    return false;

    }

    // retuns true if all but one residues in the column are gaps,
    // false otherwise

    public boolean allButOneGaps() {

	int nrOfGaps = 0;

	for (int i = 0; i < aaOcc.length; i++) {

	    if (aaOcc[i].getId() == '-')
		nrOfGaps = aaOcc[i].getOccurance();
	}

	if (nrOfGaps == columnArr.length - 1)
	    return true;

	else
	    return false;

    }

    // checks whather there is only one residue type in the column

    public boolean onlyOneResType() {

	int nrOfBig = 0;
	int indexOfBig = 0;

	for (int i = 0; i < aaOcc.length; i++) {

	    if (aaOcc[i].getOccurance() > nrOfBig) {
		nrOfBig = aaOcc[i].getOccurance();
		indexOfBig = i;
	    }
	}

	if (aaOcc[indexOfBig].getId() != '-' && nrOfBig == columnArr.length)
	    return true;

	else
	    return false;

    }

    // count the numer how many times each aa occurs in the in the column
    // the information is stored in an array 
    // of AminoAcidOccuranceChars objects that were already instantiated
    // thinking about making it a private method
    // and moving to the constructor

    private void countOccurance() {

	AminoAcidAlphabet alp = new AminoAcidAlphabet();

	char[] aAcids = alp.getAlphabet();

	for (int i = 0; i < columnArr.length; i++) {

	    for (int j = 0; j < aAcids.length; j++) {

		if (columnArr[i] == aaOcc[j].getId()) {

		    aaOcc[j].addToOccurance(1);

		    break;

		}

	    }

	}

    }

    // calculates how many residue types are  present
    // does not count gaps as residue types
    public int howManyResidueTypes() {

	int resTypeNr = 0;

	for (int i = 0; i < aaOcc.length; i++) {

	    if (aaOcc[i].getId() != '-' && aaOcc[i].getOccurance() > 0)
		resTypeNr++;
	}

	return resTypeNr;

    }

    // creates occurance objects, sets occurance to 0
    // creates 21 AminoAcidOccurance objects, the occurance is set to 0
    // id is set to one of 20 aa(one letter code) or gap(dash)

    //FIXME
    private void createOccurance() {

	AminoAcidAlphabet alp = new AminoAcidAlphabet();

	char[] aAcids = alp.getAlphabet();

	aaOcc = new AminoAcidOccurance[alp.length()];

	for (int i = 0; i < alp.length(); i++) {

	    aaOcc[i] = new AminoAcidOccurance(aAcids[i]);

	}

    }

    // gives the most common aa occurance 

    public int mostCommonOccurance() {

	int mco = 0;

	for (int i = 0; i < aaOcc.length; i++) {

	    if (aaOcc[i].getOccurance() > mco)
		mco = aaOcc[i].getOccurance();

	}

	return mco;

    }

    //calculates pairs used for Jores
    public PairCollection pairs() {

	if (this.allButOneGaps() == true) {
	    System.out.println("No pairs in this column");

	    return null;

	}

	else {

	    int diffPairs = 0;

	    int samePairs = 0;

	    int highestFreqSame = 0;

	    int highestFreqDiff = 0;

	    int highestFreq = 0;

	    for (int i = 0; i < aaOcc.length; i++) {

		if (aaOcc[i].getId() == '-' || aaOcc[i].getOccurance() == 0)
		    break;

		else {
		    samePairs++;

		    int nr1 = (aaOcc[i].getOccurance() * aaOcc[i]
			    .getOccurance())
			    - aaOcc[i].getOccurance() / 2;

		    if (nr1 > highestFreqSame)
			highestFreqSame = nr1;
		}

		for (int j = 0; j < aaOcc.length; i++) {

		    if (aaOcc[i].getId() == aaOcc[j].getId()
			    || aaOcc[j].getId() == '-')
			break;

		    else {
			diffPairs++;

			int nr2 = aaOcc[i].getOccurance()
				* aaOcc[j].getOccurance();

			if (nr2 > highestFreqDiff)
			    highestFreqDiff = nr2;
		    }

		}
	    }

	    int sumPairs = samePairs + diffPairs;

	    if (highestFreqSame > highestFreqDiff)
		highestFreq = highestFreqSame;

	    else
		highestFreq = highestFreqDiff;

	    PairCollection pair = new PairCollection(sumPairs, highestFreq);

	    return pair;

	}

    }

    // anoter way  to calculate pairs used for Jores
    // all teh special cases checked against in Jores, will be checked against here again
    // think agin about threreturns of the special cases, maybe sth better can be done
    public PairCollection pairs2() {

	PairCollection p;

	// special case #1 column is empty, only gaps

	if (this.isEmpty() == true) {
	    System.out.println("Column is empty");
	    p = new PairCollection(0, 0);
	    return p;

	}
	// special case #2	
	if (this.allButOneGaps() == true) {
	    System.out
		    .println("Only one residue in the column - no pairs formed");
	    p = new PairCollection(0, 0);
	    return p;

	}

	int pairTypes = 0;

	int mostFreqPair = 0;

	// 	special case # 3

	if (this.onlyOneResType() == true) {

	    pairTypes = 1;

	    mostFreqPair = (columnArr.length * (columnArr.length - 1)) / 2;

	}

	else {

	    pairTypes = (this.howManyResidueTypes() * (this
		    .howManyResidueTypes() - 1))
		    / 2 + this.howManyResidueTypes();

	    int big1 = 0;

	    int bigIndex = 0;

	    int big2 = 0;

	    //ArrayList<AminoAcidOccurance> arr = new ArrayList<AminoAcidOccurance>(aaOcc.length);

	    ArrayList<AminoAcidOccurance> arr = new ArrayList<AminoAcidOccurance>(
		    Arrays.asList(aaOcc));

	    assert aaOcc.length == 21;
	    //for (int i = 0; i < aaOcc.length; i++) {
	    //		arr.add(i, aaOcc[i]);
	    //	}

	    // search the list to find and remove the biggest element 

	    for (int i = 0; i < arr.size(); i++) {

		if (arr.get(i).getOccurance() > big1) {

		    big1 = arr.get(i).getOccurance();
		    bigIndex = i;
		}

	    }

	    arr.remove(bigIndex);

	    for (int i = 0; i < arr.size(); i++) {

		if (arr.get(i).getOccurance() > big2) {

		    big2 = arr.get(i).getOccurance();

		}

	    }

	    // calculates the number of time thee most frequent pair occurs
	    // if first an d second most frequent occurance equal simple product of both
	    // if not equal check pair with itself an dpairs wit second most common occurance;

	    if (big1 == big2)
		mostFreqPair = big1 * big2;

	    else {
		int most1 = big1 * big1;

		int most2 = big1 * big2;

		if (most1 > most2)
		    mostFreqPair = most1;

		else
		    mostFreqPair = most2;

	    }

	}
	p = new PairCollection(pairTypes, mostFreqPair);

	return p;

    }
}
