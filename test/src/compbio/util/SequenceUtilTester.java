package compbio.util;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class SequenceUtilTester {
	public final static String PATH= "/homes/agolicz/Eclipse/AAConservation/test/data";

	List<FastaSequence> fastlist ;

	@BeforeClass
	public void setup() throws FileNotFoundException, IOException { 
		fastlist =  SequenceUtil.readFasta(new FileInputStream(new File(PATH +"/"+ "TO1296.fasta.align")));
	}
	
	

	@Test(invocationCount=3, threadPoolSize = 3)
	public void testLoadFile() { 
		Assert.assertNotNull(fastlist); 
		Assert.assertEquals(fastlist.size(), 179);
		int length=0;
		for(FastaSequence fseq : fastlist) {
			if(length==0) {
				length= fseq.getLength();
			}
			Assert.assertEquals(length, fseq.getLength(), "Equals lenght expected");
		}

	}


	
}
