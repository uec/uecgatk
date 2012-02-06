package edu.usc.epigenome.uecgatk.YapingWalker;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import org.broadinstitute.sting.utils.GenomeLoc;

public class verboseWriter extends FormatWriter {

	public verboseWriter(File location) {
		super(location);
		// TODO Auto-generated constructor stub
	}

	public verboseWriter(OutputStream output) {
		super(output);
		// TODO Auto-generated constructor stub
	}

	public verboseWriter(File location, OutputStream output) {
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
	
	public void add(String contig, long pos, double value1, long value2) {
		// TODO Auto-generated method stub
		String valuesLine = String.format("%s\t%d\t%.2f\t%d\n",contig, pos, value1, value2);
    	try {
			mWriter.write(valuesLine);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	@Override
	public void addHeader() {
		// TODO Auto-generated method stub

	}

}
