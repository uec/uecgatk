/**
 * 
 */
package edu.usc.epigenome.uecgatk.BisSNP;

import java.util.BitSet;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import edu.usc.epigenome.uecgatk.BisSNP.BadBaseFilterBisulfite;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Mar 16, 2012 3:48:16 PM
 * 
 */
public class GATKSAMRecordFilterStorage {

	/**
	 * 
	 */

	private GATKSAMRecord GATKrecord = null;
	private boolean goodBase = false;
	private BisulfiteArgumentCollection BAC;
//	private int cytosineConvertStart = -1; //record number of cytosine converted in the reads (for minConv option in BAC)   0-based system as offset
	//private int MISMATCH_WINDOW_SIZE = 20;
	// A bitset which represents the bases of the read.  If a bit is set, then
    // the base is good; otherwise it is a bad base (as defined by the setter).
	// this is temporily use, finally we should implement indel aligner to exclude indels caused SNPs..
	//private BitSet mBitSet = null;
	
	public GATKSAMRecordFilterStorage(GATKSAMRecord GATKrecord, BisulfiteArgumentCollection BAC, ReferenceContext refContext, int offset) {
		// TODO Auto-generated constructor stub
		this.GATKrecord = GATKrecord;
		this.BAC = BAC;
		setGoodBases(offset,refContext);
	//	if(BAC.minConv > 0){
	//		setConvertStart(GATKrecord,refContext);
	//	}
		
		
		//System.err.println(cytosineConvertStart + "\t" + goodBase + "\t" + offset + "\t" + GATKrecord.getAlignmentStart() + "\t" + GATKrecord.getAlignmentEnd() + "\t" + GATKrecord.getReadNegativeStrandFlag() + "\t" + new String(GATKrecord.getReadBases()));
	}
	
	//public GATKSAMRecordFilterStorage(GATKSAMRecord GATKrecord, BisulfiteArgumentCollection BAC, ReferenceContext refContext, int offset, byte bas) {
		// TODO Auto-generated constructor stub
	//	this.GATKrecord = GATKrecord;
	//	this.BAC = BAC;
	//	setGoodBases(offset,refContext);
	//	setConvertStart(GATKrecord,refContext);
		
		//System.err.println(cytosineConvertStart + "\t" + goodBase + "\t" + offset + "\t" + GATKrecord.getAlignmentStart() + "\t" + GATKrecord.getAlignmentEnd() + "\t" + GATKrecord.getReadNegativeStrandFlag() + "\t" + new String(GATKrecord.getReadBases()) + "\t" + (char)bas);
	//}

	public GATKSAMRecordFilterStorage(GATKSAMRecord GATKrecord, BisulfiteArgumentCollection BAC, ReferenceContext refContext, int offset, boolean filterByOriginalQual) {
		// TODO Auto-generated constructor stub
		this.GATKrecord = GATKrecord;
		this.BAC = BAC;
		setGoodBases(offset,refContext, filterByOriginalQual);
	//	BitSet mismatches = BisulfiteAlignmentUtils.mismatchesInRefWindow(GATKrecord, refContext, BAC.MAX_MISMATCHES, MISMATCH_WINDOW_SIZE, BAC.sequencingMode, GATKrecord.getReadPairedFlag());
	//	if ( mBitSet == null )
	//		mBitSet = mismatches;
		
	}
	
	

    public boolean isGoodBase() {
        
    	return goodBase;
    }
    
   // public boolean isGoodConvertedBase(int offset) {
   //     if(offset >= cytosineConvertStart)
   //     	return goodBase;
   //     else
   //     	return false;
   // }
    
  //  public int getConvertStart(){
   // 	return cytosineConvertStart;
  //  }
    
	private void setGoodBases(int offset, ReferenceContext refContext) {
		byte[] quals = GATKrecord.getBaseQualities();
		
		if ( GATKrecord.getMappingQuality() >= BAC.MIN_MAPPING_QUALTY_SCORE && quals[offset] >= BAC.MIN_BASE_QUALTY_SCORE &&
	             (BAC.USE_BADLY_MATED_READS || (!BadMateFilter.hasBadMate(GATKrecord)) && !GATKrecord.getNotPrimaryAlignmentFlag()) 
	         //    && ( mBitSet == null || mBitSet.size() <= offset ? true : mBitSet.get(offset))
	             && offset >= BAC.trim5 && offset < GATKrecord.getReadLength()-BAC.trim3) {
	        	//System.out.println("bad mates");
				//if((GATKrecord.getReadPairedFlag() && GATKrecord.getProperPairFlag()))
				goodBase = true;
	     }

    }

	private void setGoodBases(int offset, ReferenceContext refContext, boolean filterByOriginalQual) {
		byte[] quals = GATKrecord.getOriginalBaseQualities();
		
		if ( GATKrecord.getMappingQuality() >= BAC.MIN_MAPPING_QUALTY_SCORE && quals[offset] >= BAC.MIN_BASE_QUALTY_SCORE &&
	             (BAC.USE_BADLY_MATED_READS || (!BadMateFilter.hasBadMate(GATKrecord)) && !GATKrecord.getNotPrimaryAlignmentFlag())
	         //    && ( mBitSet == null || mBitSet.size() <= offset ? true : mBitSet.get(offset))
	             && offset >= BAC.trim5 && offset < GATKrecord.getReadLength()-BAC.trim3) {
	        	//System.out.println("bad mates");
				//if((GATKrecord.getReadPairedFlag() && GATKrecord.getProperPairFlag()))
				goodBase = true;
	     }

    }
	/*
	private void setConvertStart(GATKSAMRecord GATKrecord, ReferenceContext refContext){
		byte[] bases = GATKrecord.getReadBases();
		if(refContext == null)
			return;
		//System.err.println(new String(refBases));
		//byte[] refBases = refContext.getBases();
		int convertedCount=0;
		if( cytosineConvertStart != -1)
			return;
		if(GATKrecord.getReadNegativeStrandFlag()){
			for(int i=0; i<bases.length;i++){
				GenomeLoc loc = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), GATKrecord.getAlignmentStart()+i);
				if( !refContext.getWindow().containsP(loc) )
					return;
				
				ReferenceContext tmpRef = new ReferenceContext(refContext.getGenomeLocParser(),loc, refContext.getWindow(),refContext.getBases());
				byte refBase = tmpRef.getBase();
				//System.err.print((char)refBase);
				if(BaseUtils.basesAreEqual(refBase, BaseUtils.G)){
					if(BaseUtils.basesAreEqual(bases[i], BaseUtils.A)){
						convertedCount++;
						if(convertedCount >= BAC.minConv){
							cytosineConvertStart=i;
							return;
						}
							
					}
				}
			}
		}
		else{
			for(int i=0; i<bases.length;i++){
				//System.err.println(GATKrecord.getAlignmentStart() + "\t" + i + "\t" + refContext.getWindow().getStart() + "\t" + refContext.getWindow().getStop());
				GenomeLoc loc = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), GATKrecord.getAlignmentStart()+i);
				if( !refContext.getWindow().containsP(loc) )
					return;
				
				ReferenceContext tmpRef = new ReferenceContext(refContext.getGenomeLocParser(),loc, refContext.getWindow(),refContext.getBases());
				byte refBase = tmpRef.getBase();
				//System.err.print((char)refBase);
				if(BaseUtils.basesAreEqual(refBase, BaseUtils.C)){
					if(BaseUtils.basesAreEqual(bases[i], BaseUtils.T)){
						convertedCount++;
						if(convertedCount >= BAC.minConv){
							cytosineConvertStart=i;
							return;
						}
							
					}
				}
			}
		}
		//System.err.println();
		
	}
	*/
	
}
