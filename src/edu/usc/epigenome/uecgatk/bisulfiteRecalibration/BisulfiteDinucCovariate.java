/**
 * 
 */
package edu.usc.epigenome.uecgatk.bisulfiteRecalibration;

import java.util.HashMap;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.walkers.recalibration.Dinuc;
import org.broadinstitute.sting.gatk.walkers.recalibration.DinucCovariate;
import org.broadinstitute.sting.gatk.walkers.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import edu.usc.epigenome.uecgatk.BisSNP.BaseUtilsMore;


/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Apr 11, 2012 8:33:47 PM
 * 
 */
public class BisulfiteDinucCovariate extends DinucCovariate {
	
	private static final byte NO_CALL = (byte) 'N';
    private static final Dinuc NO_DINUC = new Dinuc(NO_CALL, NO_CALL);
	private HashMap<Integer, Dinuc> dinucHashMap;
	private final byte[] BASES = {(byte) 'A', (byte) 'C', (byte) 'G', (byte) 'T'};
	private final byte BASE_X = (byte) 'X';
	//private final byte[] BASES_2nd_pair = {(byte) 'B', (byte) 'D', (byte) 'H', (byte) 'U', (byte) 'Y'}; // each are the base in 2nd end
	//private boolean pairedEnd=true;
	
	@Override
    public void initialize(final RecalibrationArgumentCollection RAC) {
       // BASES = {(byte) 'A', (byte) 'C', (byte) 'G', (byte) 'T', (byte) 'X'}; //X is T that converted from C (reference base is C, but read base is T)
        dinucHashMap = new HashMap<Integer, Dinuc>();
        for (byte byte1 : BASES) {
            for (byte byte2 : BASES) {
                dinucHashMap.put(Dinuc.hashBytes(byte1, byte2), new Dinuc(byte1, byte2)); // This might seem silly, but Strings are too slow
            }
            dinucHashMap.put(Dinuc.hashBytes(byte1, BASE_X), new Dinuc(byte1, BASE_X));
        }
      //  if(pairedEnd){
      //  	 for (byte byte1 : BASES_2nd_pair) {
      //           for (byte byte2 : BASES_2nd_pair) {
      //               dinucHashMap.put(Dinuc.hashBytes(byte1, byte2), new Dinuc(byte1, byte2)); // This might seem silly, but Strings are too slow
      //           }
      //       }
      //  }
       
        // Add the "no dinuc" entry too
        dinucHashMap.put(Dinuc.hashBytes(NO_CALL, NO_CALL), NO_DINUC);
    }
	
	/**
     * Takes an array of size (at least) read.getReadLength() and fills it with the covariate values for each position in the read.  TO DO: need to figure out the ref C poisiton's T is different from the other places!!!
     */
    public void getValues(final GATKSAMRecord read, final Comparable[] comparable, ReferenceContext ref) {
        final HashMap<Integer, Dinuc> dinucHashMapRef = this.dinucHashMap; //optimize access to dinucHashMap
        final int readLength = read.getReadLength();
        final boolean negativeStrand = read.getReadNegativeStrandFlag();
        boolean secondPair = false;
        if(read.getReadPairedFlag()){
        	secondPair = read.getSecondOfPairFlag();
        }
        byte[] bases = read.getReadBases();
        byte base;
        byte prevBase;
        int offset = 0;
        byte[] refBases = new byte[readLength];
		int start = read.getAlignmentStart();
				//(negativeStrand) ? read.getUnclippedEnd() : read.getAlignmentStart();
		int	end = read.getUnclippedEnd();
		GenomeLoc locBoundary = ref.getGenomeLocParser().createGenomeLoc(ref.getLocus().getContig(),start, end);
		//(negativeStrand) ? read.getAlignmentStart() : read.getUnclippedEnd();
		for(int i=0, cor=start; i<readLength; i++){
			GenomeLoc loc = ref.getGenomeLocParser().createGenomeLoc(ref.getLocus().getContig(),cor);
	        if( !ref.getWindow().overlapsP(locBoundary))
				return;
			
			ReferenceContext tmpRef = new ReferenceContext(ref.getGenomeLocParser(),loc, ref.getWindow(),ref.getBases());
			//System.err.println(tmpRef.getLocus().getStart());
			//System.err.println(tmpRef.getLocus().getStart() + "\th:" + tmpRef.getBase() + "\t" + i + "\t" + cor + "\t" + tmpRef.getWindow().getStart()+ "\t" + tmpRef.getWindow().getStop() + "\t" + ref.getWindow().getStart() + "\t" + negativeStrand + "\t" + start + "\t" + end + "\t" + new String(ref.getBases()));
			//System.err.println(tmpRef.getLocus().getStart() + "\th:" + tmpRef.getBase() + "\t" + i + "\t" + cor + "\t" + tmpRef.getWindow().getStart()+ "\t" + tmpRef.getWindow().getStop() + "\t" + ref.getWindow().getStart() + "\t" + negativeStrand + "\t" + start + "\t" + end);
			//if(tmpRef.getBase()=='C')
				//return;
			refBases[i] = tmpRef.getBase();
		//	if(negativeStrand){
				cor++;
			//}
			//else{
			//	cor++;
			//}
		}
        
        // If this is a negative strand read then we need to reverse the direction for our previous base

        if (negativeStrand) {
            bases = BaseUtils.simpleReverseComplement(bases); //this is NOT in-place
            refBases = BaseUtils.simpleReverseComplement(refBases);
        }
        comparable[0] = NO_DINUC; // No dinuc at the beginning of the read

        prevBase = bases[0];
        byte prevRefBase = refBases[0];
        byte refBase;
        offset++;
       // System.err.println(offset + "\t" + readLength + "\t" + refBases.length);
        while (offset < readLength) {
            // Note: We are using the previous base in the read, not the
            // previous base in the reference. This is done in part to be consistent with unmapped reads.
            base = bases[offset];
            refBase = refBases[offset];
            byte originalBase = bases[offset];
            if (BaseUtilsMore.isRegularBase(prevBase)) {
                if(secondPair){
                	//if(prevBase == BaseUtils.A && !negativeStrand){
                   // 	if(prevRefBase == BaseUtils.C){
                   // 		prevBase = BaseUtilsMore.X;
                   // 	}	
                   // }
                  
                    if(base == BaseUtils.A){ //
                    	if(refBase == BaseUtils.G){
                    		base = BaseUtilsMore.X;
                    	}	
                    }
                
                    //base = BaseUtilsMore.convertTo2ndEnd(base);
                    //prevBase = BaseUtilsMore.convertTo2ndEnd(prevBase);
                }
                else{
                //	if(prevBase == BaseUtils.T && !negativeStrand){
               //     	if(prevRefBase == BaseUtils.C){
               //     		prevBase = BaseUtilsMore.X;
               //     	}	
               //     }
                 //   else if(prevBase == BaseUtils.A && negativeStrand){
                 //   	if(prevRefBase == BaseUtils.G){
                 //   		prevBase = BaseUtilsMore.X;
                 //   	}
                  //  }
                    if(base == BaseUtils.T){
                    	if(refBase == BaseUtils.C){
                    		base = BaseUtilsMore.X;
                    	}	
                    }
                 //   else if(base == BaseUtils.A && negativeStrand){
                 //   	if(refBase == BaseUtils.G){
                 //   		base = BaseUtilsMore.X;
                 //   	}	
                //    }
                }
            	
            	comparable[offset] = dinucHashMapRef.get(Dinuc.hashBytes(prevBase, base));
            }
            else {
                comparable[offset] = NO_DINUC;
            }
           //System.err.println((char)prevBase + "\t" + (char)base + "\t" + (char)originalBase + "\tpreRef: " + (char)prevRefBase + "\tref: " + (char)refBase + "\t" + dinucHashMapRef.get(Dinuc.hashBytes(prevBase, base)) + "\t" + secondPair + "\t" + negativeStrand + "\t" + read.getAlignmentStart() + "\t" + offset);
            offset++;
          //  if(base == BaseUtilsMore.X && (start + offset)==7009607){
          //  	System.err.println(offset + "\t" + (char)prevBase + "\t" + (char)base + "\t" + (char)prevRefBase + "\t" + (char)refBase + "\t" + start + "\t" + end + "\t" + negativeStrand + "\t" + (start + offset) + "\n" + new String(refBases) + "\n" + new String(bases));
          //  }
            prevBase = originalBase;
            prevRefBase = refBase;
        }
        if (negativeStrand) {
            reverse(comparable);
        }
        
    }
    
 // Used to get the covariate's value from input csv file in TableRecalibrationWalker

    public Comparable getValueBisulfite(final String str) {
        byte[] bytes = str.getBytes();
        final Dinuc returnDinuc = dinucHashMap.get(Dinuc.hashBytes(bytes[0], bytes[1]));
        if (returnDinuc.compareTo(NO_DINUC) == 0) {
            return null;
        }
        return returnDinuc;
    }

    
    /**
     * Reverses the given array in place. copy from GATK as their pribate attributes
     *
     * @param array any array
     */
    private static void reverse(final Comparable[] array) {
        final int arrayLength = array.length;
        for (int l = 0, r = arrayLength - 1; l < r; l++, r--) {
            final Comparable temp = array[l];
            array[l] = array[r];
            array[r] = temp;
        }
    }
    
    
    
}
