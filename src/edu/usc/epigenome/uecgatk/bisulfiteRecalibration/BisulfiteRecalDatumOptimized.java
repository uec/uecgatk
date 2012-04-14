package edu.usc.epigenome.uecgatk.bisulfiteRecalibration;

import org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatumOptimized;
import org.broadinstitute.sting.utils.BaseUtils;

import edu.usc.epigenome.uecgatk.BisSNP.BisSNPUtils;

public class BisulfiteRecalDatumOptimized extends RecalDatumOptimized {

	public BisulfiteRecalDatumOptimized() {
		// TODO Auto-generated constructor stub
	}

	public BisulfiteRecalDatumOptimized(RecalDatumOptimized copy) {
		super(copy);
		// TODO Auto-generated constructor stub
	}

	public BisulfiteRecalDatumOptimized(long numObservations, long numMismatches) {
		super(numObservations, numMismatches);
		// TODO Auto-generated constructor stub
	}
	@Override
	public synchronized final void incrementBaseCounts( final byte curBase, final byte refBase )  {
    	
		
            //if( readBaseIndex != refBaseIndex && !(BisulfiteSnpUtil.isCytosine(refBase,false) && BisulfiteSnpUtil.isCytosine(readBase,true))) {
                //counter.novelCountsMM++;
    		if( BaseUtils.simpleBaseToBaseIndex(curBase) != BaseUtils.simpleBaseToBaseIndex(refBase) ){
    			if((BisSNPUtils.isCytosine(refBase,false) && BisSNPUtils.isCytosine(curBase,true))){
    					increment( 0, 0 );
    			}
    			else{
    				increment( 1, 1 );
    			}
            }
    		else{
    			
    			increment( 1, 0 );
    		}


        increment( 1, (BisSNPUtils.isCytosine(refBase,false) && BisSNPUtils.isCytosine(curBase,true)) ? 0 : (BaseUtils.simpleBaseToBaseIndex(curBase) == BaseUtils.simpleBaseToBaseIndex(refBase) ? 0 : 1) ); // increment takes num observations, then num mismatches
    }
}
