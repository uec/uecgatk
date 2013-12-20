package edu.usc.epigenome.uecgatk.benWalkers;

import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;

public class CpgWalkerParams {
	
	// Setting minOutputWindSize large and minScanningWindSize/maxScanningWindSize smaller
	// can produce tighter boundaries.
	public int maxScanningWindSize = 0;
	public int minOutputWindSize = 0;
	public int minScanningWindCpgs = 0;
	public boolean debug = false;
	public int minReadsForOutput = 0;
	
	
	// These only apply to variable size windows
	public boolean useVariableWindow = false;
	public int minScanningWindSize = 0;
	
	
	public MethylDbQuerier methylParams = null;

}
