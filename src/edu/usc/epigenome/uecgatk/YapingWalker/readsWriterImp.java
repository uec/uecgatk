package edu.usc.epigenome.uecgatk.YapingWalker;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import org.broadinstitute.sting.utils.GenomeLoc;

public class readsWriterImp extends FormatWriter {

	public readsWriterImp(File location) {
		super(location);
		// TODO Auto-generated constructor stub
	}

	public readsWriterImp(OutputStream output) {
		super(output);
		// TODO Auto-generated constructor stub
	}

	public readsWriterImp(File location, OutputStream output) {
		super(location, output);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void add(GenomeLoc loc, double value) {
		// TODO Auto-generated method stub

	}

	@Override
	public void add(String contig, long start, long end, double value) {
		// TODO Auto-generated method stub

	}

	@Override
	public void addHeader() {
		// TODO Auto-generated method stub

	}

	/**
	 * @param args
	 */
	
	
	public void add(String contig, long pos, byte base, int baseQ, char strand, String readID) {
		// TODO Auto-generated method stub
		String readsLine = String.format("%s\t%d\t%c\t%d\t%c\t%s\n",contig, pos, base, baseQ, strand, readID);
    	try {
			mWriter.write(readsLine);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void addCpg(String contig, long pos, byte base, int baseQ, char strand, double methyValue, String readID) {
		// TODO Auto-generated method stub
		String readsLine = String.format("%s\t%d\t%c\t%d\t%c\t%.2f\t%s\n",contig, pos, base, baseQ, strand, methyValue, readID);
    	try {
			mWriter.write(readsLine);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
