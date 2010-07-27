package compbio.util;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

import org.apache.log4j.Logger;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class SequenceUtilTester {

    public final static Logger log = Logger.getLogger(SequenceUtilTester.class);

    public final static String DATA_PATH = "test/data";

    List<FastaSequence> fastlist;

    @BeforeClass
    public void setup() throws FileNotFoundException, IOException {
	fastlist = SequenceUtil.readFasta(new FileInputStream(new File(DATA_PATH
		+ "/" + "TO1296.fasta.align")));
    }

    @Test(invocationCount = 3, threadPoolSize = 3)
    public void testLoadFile() {
	Assert.assertNotNull(fastlist);
	log.trace("Trace" + fastlist.size());
	log.debug("Debug" + fastlist.size());
	log.warn("Warn" + fastlist.size());
	log.info("Info" + fastlist.size());
	log.error("Error" + fastlist.size());
	log.fatal("Fatal" + fastlist.size());
	Assert.assertEquals(fastlist.size(), 179);
	int length = 0;
	for (FastaSequence fseq : fastlist) {
	    if (length == 0) {
		length = fseq.getLength();
	    }
	    Assert.assertEquals(length, fseq.getLength(),
		    "Equals lenght expected");
	}

    }

}
