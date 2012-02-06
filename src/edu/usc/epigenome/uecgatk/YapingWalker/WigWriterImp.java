package edu.usc.epigenome.uecgatk.YapingWalker;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import org.broadinstitute.sting.utils.GenomeLoc;

public class WigWriterImp extends FormatWriter {

	public boolean fixedSteps = false;
	private String prevChr = null;
	
	public WigWriterImp(File location) {
		super(location);
		// TODO Auto-generated constructor stub
	}

	public WigWriterImp(OutputStream output) {
		super(output);
		// TODO Auto-generated constructor stub
	}

	public WigWriterImp(File location, OutputStream output) {
		super(location, output);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void add(GenomeLoc loc, double value) {
		// TODO Auto-generated method stub
		if(prevChr == null){
			prevChr = loc.getContig();
		}
		else if(!prevChr.equalsIgnoreCase(loc.getContig())){
			if(fixedSteps){
				
			}
			else{
				addVariableStep(loc.getContig());
				prevChr = loc.getContig();
			}
		}
		String wigLine;
		if(fixedSteps){
			wigLine = String.format("%.2f\n",value);
		}
		else{
			wigLine = String.format("%d\t%.2f\n",loc.getStart(), value);
		}
		 
    	try {
			mWriter.write(wigLine);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	@Override
	public void add(String contig, long start, long end, double value) {
		// TODO Auto-generated method stub

	}

	@Override
	public void addHeader() {
		// TODO Auto-generated method stub
		try {
			mWriter.write("track type=wiggle_0 visibility=1");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private void addVariableStep(String contig) {

		this.addVariableStep(contig, 1);
	}
	
	private void addVariableStep(String contig, int span) {

		String variableStep = String.format("variableStep chrom=%s span=%d\n",contig, span);
    	try {
			mWriter.write(variableStep);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private void addFixedStep(String contig, long start, int step) {
		this.addFixedStep(contig, start, step, 1);
		
	}

	private void addFixedStep(String contig, long start, int step, int span) {

		String fixedStep = String.format("fixedStep chrom=%s start=%d step=%d span=%d\n",contig, span);
    	try {
			mWriter.write(fixedStep);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
