package edu.usc.epigenome.uecgatk.YapingWalker;

import org.broad.tribble.bed.SimpleBEDFeature;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;

public class NDRCallContext {

	private AlignmentContext context= null;
	private GenomeLoc loc = null;
	private boolean flag = false;
	private SimpleBEDFeature rodLoci = null;
	private double aveMethyInWindow = Double.NaN; //record the average GCH methylation value in the window which end in this genomic position.
	
	public NDRCallContext(AlignmentContext context, GenomeLoc loc) {
		// TODO Auto-generated constructor stub
		this.context = context;
		this.loc = loc;
	}
	
	public NDRCallContext(AlignmentContext context, GenomeLoc loc, SimpleBEDFeature rodLoci) {
		// TODO Auto-generated constructor stub
		this.context = context;
		this.loc = loc;
		this.rodLoci = rodLoci;
	}

    /** The average sequence depth inside window */
   // double aveSeqDepthWind = 0;
    
    /** The average GCH methylation value inside window */
   // double aveGchMethyWind = 0;
 
    /** The number of Gch bases called confidently (according to user threshold), either ref or other */
   // long nGchBasesCalledConfidently = 0;

    /** The sum of methylation value of Gch bases called confidently (according to user threshold), either ref or other */
  //  double sumMethyGchBasesCalledConfidently = 0;
	public AlignmentContext getRealContext(){
		return context;
	}
	
	public boolean hasRealContext(){
		if(context == null){
			return false;
		}
		return context.hasBasePileup();
	}
	
	public AlignmentContext getFakeContext(){
		
		ReadBackedPileupImpl fakePileup = new ReadBackedPileupImpl(loc);
		AlignmentContext fakeContext = new AlignmentContext(loc, fakePileup);
		return fakeContext;
	}
	
	public GenomeLoc getLoc(){
		return loc;
	}
	
	public SimpleBEDFeature getRodLoc(){
		return rodLoci;
	}
	
	public boolean getCytosinePatternFlag(){
		return flag;
	}
	
	public void setCytosinePatternFlag(boolean cytosinePatternFlag){
		this.flag = cytosinePatternFlag;
	}
	
	public void setGchMethyInWindow(double aveMethy){
		this.aveMethyInWindow = aveMethy;
	}
	
	public double getGchMethyInWindow(){
		return this.aveMethyInWindow;
	}
	
}
