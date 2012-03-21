package edu.usc.epigenome.uecgatk.BisSNP;

import java.util.BitSet;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;

import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteEnums.MethylSNPModel;


/*
 * Bis-SNP/BisSNP: It is a genotyping and methylation calling in bisulfite treated 
 * massively parallel sequencing (Bisulfite-seq and NOMe-seq) on Illumina platform
 * Copyright (C) <2011>  <Yaping Liu: lyping1986@gmail.com>

 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

public class BisulfiteAlignmentUtils extends AlignmentUtils {

	public BisulfiteAlignmentUtils() {
		// TODO Auto-generated constructor stub
	}

	/** Returns the number of mismatches in the pileup element within the given reference context in bisulfite-seq space.
    *
    * @param read          the SAMRecord
    * @param ref           the reference context
    * @param maxMismatches the maximum number of surrounding mismatches we tolerate to consider a base good
    * @param windowSize    window size (on each side) to test
    * @param sequencingMode    in Bisulfite mode, GNOMe-seq mode or Normal-seq mode
    * @param pairedend    paired-end reads or not
    * @return a bitset representing which bases are good
    */
   public static BitSet mismatchesInRefWindow(SAMRecord read, ReferenceContext ref, int maxMismatches, int windowSize, MethylSNPModel sequencingMode, boolean pairedend) {
       // first determine the positions with mismatches
       int readLength = read.getReadLength();
       BitSet mismatches = new BitSet(readLength);

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
       boolean secondPair = false;
       if(pairedend){
    	   secondPair = read.getSecondOfPairFlag();
       }

       for (int i = 0 ; i < c.numCigarElements() ; i++) {
           CigarElement ce = c.getCigarElement(i);
           int cigarElementLength = ce.getLength();
           switch ( ce.getOperator() ) {
               case M:
                   for (int j = 0; j < cigarElementLength; j++, readIndex++) {
                       // skip over unwanted bases
                       if ( currentReadPos++ < readStartPos )
                           continue;

                       // if reads extend beyond the contig end
                       if ( refIndex >= refBases.length )
                           break;

                       byte refChr = refBases[refIndex];
                       byte readChr = readBases[readIndex];
                       if ( readChr != refChr ){
                    	   if(sequencingMode == MethylSNPModel.BM || sequencingMode == MethylSNPModel.GM){ // in bisulfite conversion space, C and T will be treated as the same, but in different situation, there will be some variation..
                    		   if(!negStrand){
                        		   if(secondPair){
                        			   if(((char)refChr == 'G' && (char)readChr == 'A') || ((char)refChr == 'A' && (char)readChr == 'G')){
                                		   
                                	   }
                                	   else{
                                		   mismatches.set(readIndex);
                                	   }
                        		   }
                        		   else{
                        			   if(((char)refChr == 'C' && (char)readChr == 'T') || ((char)refChr == 'T' && (char)readChr == 'C')){
                                		   
                                	   }
                                	   else{
                                		   mismatches.set(readIndex);
                                	   }
                        		   }
                    			   
                        	   }
                        	   else{
                        		   if(secondPair){
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
           
       }

       return result;
   }
}
