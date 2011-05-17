package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import java.util.BitSet;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.NonRefDependSNPGenotypeLikelihoodsCalculationModel.MethylSNPModel;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;

public class BisulfiteAlignmentUtils extends AlignmentUtils {

	public BisulfiteAlignmentUtils() {
		// TODO Auto-generated constructor stub
	}

	/** Returns the number of mismatches in the pileup element within the given reference context in bisulfite seq space.
    *
    * @param read          the SAMRecord
    * @param ref           the reference context
    * @param maxMismatches the maximum number of surrounding mismatches we tolerate to consider a base good
    * @param windowSize    window size (on each side) to test
    * @param bisulfiteSpace    in the bisulfite conversion space
    * @return a bitset representing which bases are good
    */
   public static BitSet mismatchesInRefWindow(SAMRecord read, ReferenceContext ref, int maxMismatches, int windowSize, MethylSNPModel sequencingMode) {
       // first determine the positions with mismatches
       int readLength = read.getReadLength();
       BitSet mismatches = new BitSet(readLength);

       // it's possible we aren't starting at the beginning of a read,
       //  and we don't need to look at any of the previous context outside our window
       //  (although we do need future context)
       int readStartPos = Math.max(read.getAlignmentStart(), (int)ref.getLocus().getStart() - windowSize);
       int currentReadPos = read.getAlignmentStart();

       byte[] refBases = ref.getBases();
       int refIndex = readStartPos - (int)ref.getWindow().getStart();
       if ( refIndex < 0 ) {
           throw new IllegalStateException("When calculating mismatches, we somehow don't have enough previous reference context for read " + read.getReadName() + " at position " + ref.getLocus());
       }

       byte[] readBases = read.getReadBases();
       int readIndex = 0;

       Cigar c = read.getCigar();
       boolean negStrand = read.getReadNegativeStrandFlag();
       //System.out.println("mismatchesInRefWindow in bs");
       for (int i = 0 ; i < c.numCigarElements() ; i++) {
           CigarElement ce = c.getCigarElement(i);
           int cigarElementLength = ce.getLength();
           switch ( ce.getOperator() ) {
               case M:
                   for (int j = 0; j < cigarElementLength; j++, readIndex++) {
                       // skip over unwanted bases
                       if ( currentReadPos++ < readStartPos )
                           continue;

                       // this is possible if reads extend beyond the contig end
                       if ( refIndex >= refBases.length )
                           break;

                       byte refChr = refBases[refIndex];
                       byte readChr = readBases[readIndex];
                       if ( readChr != refChr ){
                    	   if(sequencingMode == MethylSNPModel.BM || sequencingMode == MethylSNPModel.GM){
                    		   if(!negStrand){
                        		   if(((char)refChr == 'C' && (char)readChr == 'T') || ((char)refChr == 'T' && (char)readChr == 'C')){
                            		   
                            	   }
                            	   else{
                            		   mismatches.set(readIndex);
                            	   }
                        	   }
                        	   else{
                        		   if(((char)refChr == 'G' && (char)readChr == 'A') || ((char)refChr == 'A' && (char)readChr == 'G')){
                            		   
                            	   }
                            	   else{
                            		   mismatches.set(readIndex);
                            	   }
                        	   }
                    	   }
                    	   else{
                    		   mismatches.set(readIndex);
                    	   }
                    	   
                       }
                           

                       refIndex++;
                   }
                   break;
               case I:
               case S:
                   readIndex += cigarElementLength;
                   break;
               case D:
               case N:
                   if ( currentReadPos >= readStartPos )
                       refIndex += cigarElementLength;
                   currentReadPos += cigarElementLength;
                   break;
               case H:
               case P:
                   break;
           }
       }

       // all bits are set to false by default
       BitSet result = new BitSet(readLength);

       int currentPos = 0, leftPos = 0, rightPos;
       int mismatchCount = 0;

       // calculate how many mismatches exist in the windows to the left/right
       for ( rightPos = 1; rightPos <= windowSize && rightPos < readLength; rightPos++) {
           if ( mismatches.get(rightPos) )
               mismatchCount++;
       }
       if ( mismatchCount <= maxMismatches )
           result.set(currentPos);
       else{
    	   //System.out.println("1\t" + currentPos);
       }
       //System.out.println("3\t" + mismatchCount);
       // now, traverse over the read positions
       while ( currentPos < readLength ) {
           // add a new rightmost position
           if ( rightPos < readLength && mismatches.get(rightPos++) )
               mismatchCount++;
           // re-penalize the previous position
           if ( mismatches.get(currentPos++) )
               mismatchCount++;
           // don't penalize the current position
           if ( mismatches.get(currentPos) )
               mismatchCount--;
           // subtract the leftmost position
           if ( leftPos < currentPos - windowSize && mismatches.get(leftPos++) )
                   mismatchCount--;

           if ( mismatchCount <= maxMismatches )
               result.set(currentPos);
           else{
        	   //System.out.println("2\t" + currentPos);
           }
       }

       //System.out.println("4\t" + mismatchCount);
       return result;
   }
}
