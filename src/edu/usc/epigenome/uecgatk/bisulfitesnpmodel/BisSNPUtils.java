package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeSet;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;

import edu.usc.epigenome.uecgatk.YapingWalker.NDRargumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisulfiteGenotyperEngine.BadBaseFilterBisulfite;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisulfiteGenotyperEngine.OUTPUT_MODE;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisulfiteSNPGenotypeLikelihoodsCalculationModel.methyStatus;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.NonRefDependSNPGenotypeLikelihoodsCalculationModel.MethylSNPModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecordFilter;

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

/*
 * provide easy access for some function inside BisSNP.
 *  e.g. judge if it is a type of cytosine pattern(HCG-2, or GCH-1) from the given pileup and corresponding reference seq, dbSNP. 
 * should provide methylation value, otherwise will use flat methylation value; and also likelihood ratio criteria, bisulfite conversion rate
 */

public class BisSNPUtils {
	
	private double FLAT_METHY_STATUS = 0.5;
	private NDRargumentCollection NAC;
	private BisulfiteDiploidSNPGenotypePriors genotypePriors;
	private BisulfiteArgumentCollection BAC;
	
	public BisSNPUtils(BisulfiteArgumentCollection BAC){
		this.NAC = new NDRargumentCollection();
		this.genotypePriors = new BisulfiteDiploidSNPGenotypePriors();
	}
	
	public BisSNPUtils(NDRargumentCollection NAC){
		this.NAC = NAC;
		this.genotypePriors = new BisulfiteDiploidSNPGenotypePriors();
	}
	
	public boolean isGch(ReadBackedPileup pileup, RefMetaDataTracker tracker,ReferenceContext ref, double methyStatus){
		return checkCytosineStatus("GCH-2", pileup, tracker, ref, genotypePriors, BAC, methyStatus);
		
	}
	
	public boolean isHcg(ReadBackedPileup pileup, RefMetaDataTracker tracker,ReferenceContext ref, BisulfiteDiploidSNPGenotypePriors priors, 
			BisulfiteArgumentCollection bac, double methyStatus){
		return checkCytosineStatus("HCG-2", pileup, tracker, ref, genotypePriors, BAC, methyStatus);
		
	}
	
	public boolean isWcg(ReadBackedPileup pileup, RefMetaDataTracker tracker,ReferenceContext ref, BisulfiteDiploidSNPGenotypePriors priors, 
			BisulfiteArgumentCollection bac, double methyStatus){
		return checkCytosineStatus("WCG-2", pileup, tracker, ref, genotypePriors, BAC, methyStatus);
		
	}
	
	public boolean isCpg(ReadBackedPileup pileup, RefMetaDataTracker tracker,ReferenceContext ref, BisulfiteDiploidSNPGenotypePriors priors, 
			BisulfiteArgumentCollection bac, double methyStatus){
		return checkCytosineStatus("CG-1", pileup, tracker, ref, genotypePriors, BAC, methyStatus);
		
	}
	
	public boolean isCph(ReadBackedPileup pileup, RefMetaDataTracker tracker,ReferenceContext ref, BisulfiteDiploidSNPGenotypePriors priors, 
			BisulfiteArgumentCollection bac, double methyStatus){
		return checkCytosineStatus("CH-1", pileup, tracker, ref, genotypePriors, BAC, methyStatus);
		
	}
	
	public boolean isCytosine(ReadBackedPileup pileup, RefMetaDataTracker tracker,ReferenceContext ref, BisulfiteDiploidSNPGenotypePriors priors, 
			BisulfiteArgumentCollection bac, double methyStatus){
		return checkCytosineStatus("C-1", pileup, tracker, ref, genotypePriors, BAC, methyStatus);
		
	}
	
	public boolean isCytosineType(ReadBackedPileup pileup, RefMetaDataTracker tracker,ReferenceContext ref, BisulfiteDiploidSNPGenotypePriors priors, 
			BisulfiteArgumentCollection bac, double methyStatus, String cytosineTypeToCheck){ //should be "CH-1" style..
		return checkCytosineStatus(cytosineTypeToCheck, pileup, tracker, ref, genotypePriors, BAC, methyStatus);
		
	}
	
	public boolean checkCytosineStatus(String cytosinePattern, ReadBackedPileup pileup, RefMetaDataTracker tracker,ReferenceContext ref, BisulfiteDiploidSNPGenotypePriors priors, 
			BisulfiteArgumentCollection bac, String methySummaryFileName){
		
		FileReader methySummaryFile=null;
		try {
			methySummaryFile = new FileReader(methySummaryFileName);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		BufferedReader br = new BufferedReader(methySummaryFile);
		String line;
		double methyStatus = 0.5;
		try {
			while( (line = br.readLine()) != null){
				if(line.startsWith(cytosinePattern)){
					methyStatus = Double.parseDouble(line.split("\t")[0].split(":")[1]);
					break;
				}
					
				
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return checkCytosineStatus(cytosinePattern, pileup, tracker, ref, priors, bac, methyStatus); 
	}
	
	public boolean checkCytosineStatus(String cytosinePattern, ReadBackedPileup pileup, RefMetaDataTracker tracker,ReferenceContext ref, BisulfiteDiploidSNPGenotypePriors priors, 
			BisulfiteArgumentCollection bac, double methyStatus){
		//double[] cytosineMethyStatus = new double[2]; // 0: methy status in positive strand; 1: methy status in negative strand;
		String[] tmpKey = cytosinePattern.split("-");
		int maxCytosineLength = Math.max(tmpKey[0].length()-Integer.parseInt(tmpKey[1]),Integer.parseInt(tmpKey[1]));
		//BisulfiteDiploidSNPGenotypeLikelihoods maxGL = null;
		GenomeLoc location = pileup.getLocation();
		String contig = location.getContig();
		int position = location.getStart();
		double[] tmpMethy = new double[2];
		tmpMethy[0] = FLAT_METHY_STATUS;
		tmpMethy[1] = FLAT_METHY_STATUS;
		HashMap<Integer,methyStatus> cytosineAndAdjacent = new HashMap<Integer,methyStatus>();
		//check adjacent position likelihood
		for(int i = 0 - maxCytosineLength; i <= maxCytosineLength; i++){
			GenomeLoc loc = ref.getGenomeLocParser().createGenomeLoc(contig, position + i );
			if(i == 0)
				continue;
			List<SAMRecord> reads =  new ArrayList<SAMRecord>();;
			List<Integer> elementOffsets = new ArrayList<Integer>();

			for ( PileupElement p : pileup ) {
					int elementOffset = i + p.getOffset();
					if(elementOffset < 0 || elementOffset > p.getRead().getReadLength()-1)
						continue;
					elementOffsets.add(elementOffset);
					reads.add(p.getRead());
			}
			ReadBackedPileup tmpPileup = new ReadBackedPileupImpl(loc,reads,elementOffsets);
			
		//	ReferenceContext tmpRef = new ReferenceContext(ref.getGenomeLocParser(),loc, ref.getWindow(),ref.getBases());
			 
			ReferenceContext tmpRef = new ReferenceContext(ref.getGenomeLocParser(),loc, ref.getWindow(), ref.getBases());
			BisulfiteDiploidSNPGenotypeLikelihoods tmpGL = new BisulfiteDiploidSNPGenotypeLikelihoods(tracker, tmpRef, (BisulfiteDiploidSNPGenotypePriors)priors, bac, tmpMethy.clone());
			
				
			tmpGL.setPriors(tracker, tmpRef, bac.heterozygosity, bac.novelDbsnpHet, bac.validateDbsnpHet, loc);
			if(position == bac.testLocus){
            	System.err.println("i: " + i + "\ttmpRef: " + tmpRef.getBase());
            	tmpGL.VERBOSE = true;
            }
			int nGoodBases = tmpGL.add(tmpPileup, true, true);
            if ( nGoodBases == 0 )
                continue;
            double[] posteriorNormalized = normalization(tmpGL.getPosteriors(),tmpGL.getLikelihoods());
            
            getBestGenotypeFromPosterior(posteriorNormalized, cytosineAndAdjacent, i ,position);
            
		}
		BisulfiteDiploidSNPGenotypeLikelihoods tmpGL = new BisulfiteDiploidSNPGenotypeLikelihoods(tracker, ref, (BisulfiteDiploidSNPGenotypePriors)priors, bac, tmpMethy.clone());
		tmpGL.setPriors(tracker, ref, bac.heterozygosity, bac.novelDbsnpHet, bac.validateDbsnpHet, location);
		
			//Double[] value = cts.cytosineListMap.get(cytosineType);
			tmpMethy[0] = methyStatus;
			tmpMethy[1] = methyStatus;
			int cytosinePos = Integer.parseInt(tmpKey[1]);
			
			
          //  double adjacentCytosineSeqLikelihood = 0;
			//double adjacentCytosineSeqLikelihoodReverseStrand = 0;
			int i = 1;
			int countMatchedOnFwd = 0;
			int countMatchedOnRvd = 0;
            for(byte base : tmpKey[0].getBytes()){
            	int pos = i - cytosinePos;
            	i++;
            	if(pos == 0)
            		continue;
            	methyStatus tmpMethyStatus = cytosineAndAdjacent.get(pos);
            	if(tmpMethyStatus == null){
            		break;
            	}
            	else if(tmpMethyStatus.genotype == null){
	            	break;
	            }
	            else {
	            	 
	                 if(tmpMethyStatus.ratio < bac.STANDARD_CONFIDENCE_FOR_CALLING){
	                 		break;
	                 }

	            	 if(tmpMethyStatus.genotype.isHet()){
	 	            	break;
	 	             }	
	 	             else{
	 	            	 
	 	            	if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, tmpMethyStatus.genotype.base1)){
	 	            		countMatchedOnFwd++;
	 	            	//	adjacentCytosineSeqLikelihood += tmpMethyStatus.ratio;
	 	            	}
	 	            	
	 	             }
	            }
            	
	            if(position == bac.testLocus){
	            	System.err.println("base: " + (char)base + "\tgenotype: " + (char)tmpMethyStatus.genotype.base1 + "\tcytosinePos: " + cytosinePos + "\tratio: " + tmpMethyStatus.ratio);
	            }
	            
            }
            i = 1;
            for(byte base : tmpKey[0].getBytes()){
            	int pos = cytosinePos - i;
            	i++;
            	if(pos == 0)
            		continue;
            	methyStatus tmpMethyStatus = cytosineAndAdjacent.get(pos);
            	if(tmpMethyStatus == null){
            		break;
            	}
            	else if(tmpMethyStatus.genotype == null){
	            	break;
	            }
            	else {
	            	 
	                 if(tmpMethyStatus.ratio < bac.STANDARD_CONFIDENCE_FOR_CALLING){
	                 	break;
	                 }
	                 	
	            	 if(tmpMethyStatus.genotype.isHet()){
	 	            	break;
	 	             }	
	 	             else{
	 	            	 
	 	            	if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(base, BaseUtilsMore.iupacCodeComplement(tmpMethyStatus.genotype.base1))){
		            		countMatchedOnRvd++;
		            	//	adjacentCytosineSeqLikelihoodReverseStrand += tmpMethyStatus.ratio;
		            	}
	 	            	
	 	             }
	            }
	            
	            if(position == bac.testLocus){
	            	System.err.println("base: " + (char)base + "\tgenotype: " + (char)BaseUtilsMore.iupacCodeComplement(tmpMethyStatus.genotype.base1) + "\treveser: " + (char)base + "\tcytosinePos: " + cytosinePos + "\tratio: " + tmpMethyStatus.ratio );
	            }
	            
            }
            if((countMatchedOnFwd < tmpKey[0].length() - 1) && (countMatchedOnRvd < tmpKey[0].length() - 1))
            	return false;
            
            
            	tmpGL.clearLikelihoods(tmpMethy.clone());
            	if(position == bac.testLocus){
            		tmpGL.VERBOSE = true;
            		System.err.println("cytosineType: " + cytosinePattern );
            	}
            	
     			int nGoodBases = tmpGL.add(pileup, true, true);
                 if ( nGoodBases == 0 )
                	 return false;
                 double[] posteriorNormalized = normalization(tmpGL.getPosteriors(),tmpGL.getLikelihoods());
                 
                 getBestGenotypeFromPosterior(posteriorNormalized, cytosineAndAdjacent, 0 ,position);

            
           
            
            methyStatus tmpMethyStatus = cytosineAndAdjacent.get(0);
        	if(tmpMethyStatus == null){
        		return false;
        	}
        	else if(tmpMethyStatus.genotype == null){
        		return false;
            }
            else if(tmpMethyStatus.genotype.isHet()){
            	return false;
            		
            }	
            else{
            	
            	if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(BaseUtils.C, tmpMethyStatus.genotype.base1)){
            		
                     	if(tmpMethyStatus.ratio < bac.STANDARD_CONFIDENCE_FOR_CALLING){
                     		return false;
                     	}
                    
            		countMatchedOnFwd++;
            	//	adjacentCytosineSeqLikelihood += tmpMethyStatus.ratio;
            	}
            	else if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(BaseUtils.C, BaseUtilsMore.iupacCodeComplement(tmpMethyStatus.genotype.base1))){
            		
                     	if(tmpMethyStatus.ratio < bac.STANDARD_CONFIDENCE_FOR_CALLING){
                     		return false;
                     	}
                    
            		countMatchedOnRvd++;
            	//	adjacentCytosineSeqLikelihoodReverseStrand += tmpMethyStatus.ratio;
            	}
            	else{
            		return false;
                	
            	}
            	
            }
            
            if(countMatchedOnFwd >= tmpKey[0].length()){
				return true;

			}
			else if(countMatchedOnRvd >= tmpKey[0].length()){
				return true;
				
					
			}
			if(position == bac.testLocus){
            	System.err.println("countMatchedOnFwd: " + countMatchedOnFwd + "\tcountMatchedOnRvd: " + countMatchedOnRvd);
            } 
		
			return false;
		
	}
	
	public void getBestGenotypeFromPosterior(double[] posterior,HashMap<Integer,methyStatus> cytosineAdjacent, int key, int location){
		double maxCount = Double.NEGATIVE_INFINITY;
        double secondMaxCount = Double.NEGATIVE_INFINITY;
        methyStatus tmpMethyStatus = new methyStatus();
        tmpMethyStatus.genotype = null;
        tmpMethyStatus.ratio = 0.0;
        DiploidGenotype bestGenotype = DiploidGenotype.createHomGenotype(BaseUtils.A);
 
        for ( DiploidGenotype g : DiploidGenotype.values() ){
			if(posterior[g.ordinal()] > maxCount){
				secondMaxCount = maxCount;
				maxCount = posterior[g.ordinal()];
				
				bestGenotype = g;
			}
			else if (posterior[g.ordinal()] > secondMaxCount && posterior[g.ordinal()] <= maxCount){
	            	secondMaxCount = posterior[g.ordinal()];
	        }
		}
        tmpMethyStatus.ratio = 10 * (maxCount - secondMaxCount);
      //  if(location == bac.testLocus){
      //  	System.err.println("maxCount: " + maxCount + "\tsecondMaxCount: " + secondMaxCount + "\tratio: " + tmpMethyStatus.ratio + "\tgenotype: " + bestGenotype);
      //  	for(double poster : posterior){
      //  		System.err.println(poster);
      //  	}
      //  }
        tmpMethyStatus.genotype = bestGenotype;
        cytosineAdjacent.put(key, tmpMethyStatus);
        
      
	}
	
	//normalize of posterior
	public double[] normalization(double[] logPosterior, double[] logLikilyhood){
		double sum = 0;
		double[] returnLikilyhood = logPosterior.clone();
		for(int i = 0; i < logLikilyhood.length; i++){
			sum += Math.pow(10,logLikilyhood[i]);
		}
		sum = Math.log10(sum);
		for(int j = 0; j < logLikilyhood.length; j++){
			returnLikilyhood[j] = returnLikilyhood[j] - sum;
		}
		return returnLikilyhood;
	}
	
	public AlignmentContext getFilteredAndStratifiedContexts(BisulfiteArgumentCollection BAC, ReferenceContext refContext, AlignmentContext rawContext) {
		BadBaseFilterBisulfite badReadPileupFilter = new BadBaseFilterBisulfite(refContext, BAC);

		//AlignmentContext stratifiedContexts = null;
       
        if ( !rawContext.hasExtendedEventPileup() ) {

            byte ref = refContext.getBase();
            if ( !BaseUtils.isRegularBase(ref) )
                return null;

            if ( !filterPileupBisulfite(rawContext, badReadPileupFilter, BAC) )
                return null;
        }

        return rawContext;
    }

	protected boolean filterPileupBisulfite(AlignmentContext stratifiedContexts, BadBaseFilterBisulfite badBaseFilter, BisulfiteArgumentCollection BAC) {
	        	int numDeletions = 0, pileupSize = 0;

	        	if(!stratifiedContexts.hasBasePileup())
	        		return false;
	       
	            ReadBackedPileup pileup = stratifiedContexts.getBasePileup();
	            for ( PileupElement p : pileup ) {
	                final SAMRecord read = p.getRead();
	                if(read.getDuplicateReadFlag()){ //get rid of duplicate reads
	                	continue;
	                }
	                if ( p.isDeletion() ) {
	                    // if it's a good read, count it
	                    if ( read.getMappingQuality() >= BAC.MIN_MAPPING_QUALTY_SCORE &&
	                         (BAC.USE_BADLY_MATED_READS || !BadMateFilter.hasBadMate(read)) )
	                        numDeletions++;
	                } else {
	                    if ( !(read instanceof GATKSAMRecord) )
	                        throw new ReviewedStingException("The BisulfiteGenotyper currently expects GATKSAMRecords, but instead saw a " + read.getClass());
	                    GATKSAMRecord GATKrecord = (GATKSAMRecord)read;
	                    GATKrecord.setGoodBases(badBaseFilter, true);
	                    if ( GATKrecord.isGoodBase(p.getOffset()) )
	                        pileupSize++;
	                }
	            }
	        

	        if ( BAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES )
	            return true;

	        // if no coverage?
	        if ( pileupSize == 0 )
	            return false;

	        // too many deletions in the pileup?
	        if ( (BAC.MAX_DELETION_FRACTION >=0 && BAC.MAX_DELETION_FRACTION <=1.0 ) &&
	                (double)numDeletions / (double)(pileupSize + numDeletions) > BAC.MAX_DELETION_FRACTION )
	            return false;

	        return true;
	    }

	 //copy from GATK, since it is not public class there
	protected class BadBaseFilterBisulfite implements GATKSAMRecordFilter {
        private ReferenceContext refContext;
        private final BisulfiteArgumentCollection BAC;
        private final int MISMATCH_WINDOW_SIZE = 20;

        public BadBaseFilterBisulfite(ReferenceContext refContext, BisulfiteArgumentCollection BAC) {
            this.refContext = refContext;
            this.BAC = BAC;
        }

        @Override
        public BitSet getGoodBases(final GATKSAMRecord record) {
            BitSet bitset = new BitSet(record.getReadLength());

            // if the mapping quality is too low or the mate is bad, we can just zero out the whole read and continue
            if ( record.getMappingQuality() < BAC.MIN_MAPPING_QUALTY_SCORE ||
                 (!BAC.USE_BADLY_MATED_READS && BadMateFilter.hasBadMate(record)) ) {
            	//System.out.println("bad mates");
            	return bitset;
            }
            
            byte[] quals = record.getBaseQualities();
            for (int i = 0; i < quals.length; i++) {
                if ( quals[i] >= BAC.MIN_BASE_QUALTY_SCORE )
                    bitset.set(i);
            }

            // if a read is too long for the reference context, extend the context (being sure not to extend past the end of the chromosome)
            if ( record.getAlignmentEnd() > refContext.getWindow().getStop() ) {
                GenomeLoc window = refContext.getWindow();
                byte[] bases = refContext.getBases();
                StringUtil.toUpperCase(bases);
                refContext = new ReferenceContext(refContext.getGenomeLocParser(),refContext.getLocus(), window, bases);
            }            

            BitSet mismatches;
            boolean paired = record.getReadPairedFlag();
            if(BAC.sequencingMode == MethylSNPModel.BM || BAC.sequencingMode == MethylSNPModel.GM){
            	mismatches = BisulfiteAlignmentUtils.mismatchesInRefWindow(record, refContext, BAC.MAX_MISMATCHES, MISMATCH_WINDOW_SIZE, BAC.sequencingMode, paired);
            }
            else{
            	mismatches = AlignmentUtils.mismatchesInRefWindow(record, refContext, BAC.MAX_MISMATCHES, MISMATCH_WINDOW_SIZE);
            }
            if ( mismatches != null )
                bitset.and(mismatches);

            return bitset;
        }
    }
	
	public static boolean usableBase(PileupElement p, boolean ignoreBadBases) {
        // ignore deletions, Q0 bases, and filtered bases
        if ( p.isDeletion() ||
                p.getQual() == 0 ||
                (p.getRead() instanceof GATKSAMRecord &&
                 !((GATKSAMRecord)p.getRead()).isGoodBase(p.getOffset())) )
            return false;

        return ( !ignoreBadBases || !badBase(p.getBase()) );
    }
	
	public static boolean badBase(byte observedBase) {
        return BaseUtils.simpleBaseToBaseIndex(observedBase) == -1;
    }
	
	public static ReadBackedPileup getDownsampledPileup(ReadBackedPileup pileup, int desiredCov){
		if ( pileup.size() <= desiredCov )
            return pileup;

        // randomly choose numbers corresponding to positions in the reads list
        Random generator = new Random();
        TreeSet<Integer> positions = new TreeSet<Integer>();
        for ( int i = 0; i < desiredCov; /* no update */ ) {
            if ( positions.add(generator.nextInt(pileup.size())) )
                i++;
        }
		GenomeLoc loc = pileup.getLocation();
		List<SAMRecord> reads =  new ArrayList<SAMRecord>();;
		List<Integer> elementOffsets = new ArrayList<Integer>();
		int i = 0;
		for ( PileupElement p : pileup ) {
			if(positions.contains(i)){
				int elementOffset = p.getOffset();
				if(elementOffset < 0 || elementOffset > p.getRead().getReadLength()-1)
					continue;
				elementOffsets.add(elementOffset);
				reads.add(p.getRead());
			}
				i++;
		}
		ReadBackedPileup downsampledPileup = new ReadBackedPileupImpl(loc,reads,elementOffsets);
		return downsampledPileup;
	}
	
	public static String getRefGenomeVersion(){
		return BisulfiteGenotyper.getBAC().referenceGenomeVer;
	}
	
}
