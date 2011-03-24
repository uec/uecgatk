package edu.usc.epigenome.uecgatk;

import java.io.File;
import java.io.OutputStream;

import org.broadinstitute.sting.utils.wiggle.WiggleWriter;

public class WiggleWriterReducible extends WiggleWriter {

	protected File fFile = null;
	
	public WiggleWriterReducible(File arg0) {
		super(arg0);
		fFile = arg0;
	}

	public WiggleWriterReducible(OutputStream out) throws Exception {
		super(out);
		
		throw new Exception("WiggleWriterReducible must be called with a Filename");
	}
	
	public static WiggleWriterReducible merge(WiggleWriterReducible a, WiggleWriterReducible b, File newFile)
	{
		return a;
	}

}
