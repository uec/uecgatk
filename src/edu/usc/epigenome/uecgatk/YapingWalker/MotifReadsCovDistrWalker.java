/**
 * Give the motif frequency, read coverage average and standard variation summary. 
 * Also output read coverage file under this motif. Only support single motif yet.
 * // 1: reference genome size (exclude 'N/n' position); 2: motif freq in reference genome; 3: motif coverage in reference genome; 
            //4: motif freq in CGI(or custom interval); 5: motif coverage in CGI(or custom interval); 
            //6: motif freq not in CGI(or custom interval); 7: motif coverage not in CGI(or custom interval); 8: motif methylation level
            // 9: motif num C in total; 10: motif num T in total; 11: strand bias(CT strand_num/ GA strand_num); 12: 1st end/2nd end bias
            //interpret the following later..
            //13: mismatch rate in total; 14: mismatch rate in 1st end; 15: mismatch rate in 2nd end. 
 */
package edu.usc.epigenome.uecgatk.YapingWalker;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import net.sf.picard.reference.IndexedFastaSequenceFile;

import org.apache.commons.math3.stat.descriptive.SynchronizedMultivariateSummaryStatistics;
import org.broad.tribble.bed.BEDFeature;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Downsample;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import edu.usc.epigenome.uecgatk.BisSNP.BaseUtilsMore;
import edu.usc.epigenome.uecgatk.YapingWalker.MotifFreqInGenomeWalker.Datum;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Dec 10, 2012 10:31:53 AM
 * 
 */
@Reference(window=@Window(start=-500,stop=500))
@Requires({DataSource.REFERENCE, DataSource.REFERENCE_BASES, DataSource.READS})
@By(DataSource.REFERENCE)
@Downsample(by=DownsampleType.NONE)
public class MotifReadsCovDistrWalker extends LocusWalker<Long, Long> implements TreeReducible<Long> {

	@Input(fullName = "motif_to_search", shortName = "motif", doc = "motif pattern to search in provided genome (only support single motif yet)", required = true)
	public String motif = null;
	
	@Input(fullName = "cgi_file", shortName = "cgi", doc = "Give the CGI bed file for the coverage estimation", required = false)
	public RodBinding<BEDFeature> cgi = null;
	
	@Input(fullName = "dbsnp", shortName = "D", doc = "dbSNP file", required = false)
	public RodBinding<VariantContext> dbsnp;
	
	@Output(fullName = "motif_read_coverage", shortName = "covDistr", doc = "coverage detail under the provided motif, could be used to plot histgram further", required = true) //2nd column is whether or not in CGI, 1: yes; 0: no.
	public PrintStream covDistr = null;
	
	
	private SynchronizedMultivariateSummaryStatistics summary;

	/**
	 * 
	 */
	public void initialize(){
	
            summary = new SynchronizedMultivariateSummaryStatistics(14, true); 
            // 1: reference genome size (exclude 'N/n' position); 2: motif number in reference genome; 3: motif coverage in reference genome; 
            //4: reference genome size in CGI(or custom interval); 5: motif num in CGI(or custom interval); 6: motif coverage in CGI(or custom interval); 
            //7: motif coverage not in CGI(or custom interval);
            //8: strand bias(motif strand_num/ not_motif strand_num); 9: strand bias(positive strand_num/ negative strand_num); 
            // 10. motif number that has coverage..
            //interpret the following later..
            //11: 1st end/2nd end bias(1st end/total); 12: mismatch rate sum in 1st end; 13: mismatch rate sum in 2nd end; 14: loci has coverage and not in dbSNP position
   
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.LocusWalker#map(org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker, org.broadinstitute.sting.gatk.contexts.ReferenceContext, org.broadinstitute.sting.gatk.contexts.AlignmentContext)
	 */

	public Long map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
		double[] value = new double[14];
		for(int i=0; i<value.length;i++){
			value[i] = 0;
		}

		boolean inCgi = false;
		if(BaseUtils.isRegularBase(ref.getBase())){
			value[0]=1;
			inCgi = isInCgi(tracker);
			if(inCgi){
				value[3]=1;
			}
			else{
				value[3]=0;
			}
		}	
		else{
			return null;
		}
		
		if (checkPattern(motif, ref, false) ){
			
			value[1]=1;
			value[2]=context.size();
			if(inCgi){
				value[4]=1;
				value[5]=context.size();
			}
			else{
				value[4]=0;
				value[6]=context.size();
			}

			if(context.hasBasePileup() && !context.getBasePileup().isEmpty()){
				value[9] = 1;
				double numPos = context.getBasePileup().getPositiveStrandPileup().depthOfCoverage();
				double numNeg = context.getBasePileup().getNegativeStrandPileup().depthOfCoverage();
				value[7] = numPos/(numNeg+numPos);
				value[8] = numPos/(numNeg+numPos);
				value[10] = firstEndFraction(context);
				if(tracker.getValues(dbsnp).isEmpty()){
					value[13] = 1;
					value[11] = Double.compare(value[10], 0.0) == 0 ? 0 : mismatchFraction(ref, context, true) / (value[10] * value[2]);
					value[12] = Double.compare(value[10], 1.0) == 0 ? 0 : mismatchFraction(ref, context, false) / ((1-value[10]) * value[2]);
				}

				System.err.println(value[8] + "\t" + numPos + "\t" + numNeg);
			}
			covDistr.println(value[2] + "\t" + value[4]);
		}
		else if(checkPattern(motif, ref, true)){

			value[1]=1;
			value[2]=context.size();
			if(inCgi){
				value[4]=1;
				value[5]=context.size();
			}
			else{
				value[4]=0;
				value[6]=context.size();
			}

			if(context.hasBasePileup() && !context.getBasePileup().isEmpty()){
				value[9] = 1;
				double numPos = context.getBasePileup().getPositiveStrandPileup().depthOfCoverage();
				double numNeg = context.getBasePileup().getNegativeStrandPileup().depthOfCoverage();
				value[7] = numNeg/(numNeg+numPos);
				value[8] = numPos/(numNeg+numPos);
				value[10] = firstEndFraction(context);
				System.err.println(value[7] + "\t" + numPos + "\t" + numNeg);
				if(tracker.getValues(dbsnp).isEmpty()){
					value[13] = 1;
					value[11] = Double.compare(value[10], 0.0) == 0 ? 0 : mismatchFraction(ref, context, true) / (value[10] * value[2]);
					value[12] = Double.compare(value[10], 1.0) == 0 ? 0 : mismatchFraction(ref, context, false) / ((1-value[10]) * value[2]);
				}
			}
			covDistr.println(value[2] + "\t" + value[4]);
		}
		summary.addValue(value);
		
		 return null;
	}
	
	public void onTraversalDone(Long result) {
		double[] sum = summary.getSum();

		double[] sd = summary.getStandardDeviation();
		logger.info("Genome size(Only count regular bases): " + (int)sum[0]);
		logger.info("Genome size(Only count regular bases, in Given Interval(e.g. CGI)): " + (int)sum[3]);
		logger.info("Motif Pattern " + motif + " number: " + (int)sum[1]);
		logger.info("Motif frequency (%): " + String.format("%.3f", 100*sum[1]/sum[0]));
		logger.info("Motif coverage mean: " +  String.format("%.1f", sum[2]/sum[1]));
		logger.info("Motif coverage standard deviation: " + sd[2]);
		
		logger.info("Motif Pattern " + motif + " number (in CGI): " + (int)sum[4]);
		logger.info("Motif frequency (%) (in CGI) : " + String.format("%.3f", 100*sum[4]/sum[3]));
		logger.info("Motif coverage mean (in CGI) : " + String.format("%.1f", sum[5]/sum[4]));
		logger.info("Motif coverage standard deviation (in CGI) : " + sd[5]);
		
		logger.info("Motif Pattern " + motif + " number (NOT in CGI): " + (int)(sum[1] - sum[4]));
		logger.info("Motif frequency (%) (NOT in CGI) : " + String.format("%.3f", 100*(sum[1] - sum[4])/(sum[0] - sum[3])));
		logger.info("Motif coverage mean (NOT in CGI) : " + String.format("%.1f", sum[6]/(sum[1] - sum[4])));
		logger.info("Motif coverage standard deviation (NOT in CGI) : " + sd[6]);
		
		logger.info("Motif coverage ratio (in CGI vs. NOT in CGI) : " + String.format("%.3f", (sum[2]/sum[1])/(sum[6]/(sum[1] - sum[4]))));
		
		logger.info("Motif strand bias mean (Pos/Neg strand) : " + String.format("%.3f", sum[8]/sum[9]));
		logger.info("Motif strand bias standard deviation (Pos/Neg strand) : " + sd[8]);
		
		logger.info("Motif strand bias mean (Motif_strand/No_Motif_strand) : " + String.format("%.3f", sum[7]/sum[9]));
		logger.info("Motif strand bias standard deviation (Motif_strand/No_Motif_strand) : " + sd[7]);
		
		logger.info("Motif 1st end ratio mean (1st_end/total_reads) : " + String.format("%.3f", sum[10]/sum[9]));
		logger.info("Motif 1st end ratio standard deviation (Pos/Neg strand) : " + sd[10]);

		logger.info("Motif 1st end mismatch rate mean (%) : " + String.format("%.3f", 100*sum[11]/sum[13]));
		logger.info("Motif 1st end ratio standard deviation (Pos/Neg strand) : " + sd[11]);
		
		logger.info("Motif 2nd end mismatch rate mean (%) : " + String.format("%.3f", 100*sum[12]/sum[13]));
		logger.info("Motif 2nd end ratio standard deviation (Pos/Neg strand) : " + sd[12]);
		
		logger.info("Finished!");
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduceInit()
	 */

	public Long reduceInit() {
		// TODO Auto-generated method stub
		return null;
	}

	
	private boolean checkPattern(String motif, ReferenceContext ref, boolean negStrand){
		byte[] refBytes = new byte[motif.length()];
		byte[] motifSeq = motif.getBytes();
		if(negStrand)
			motifSeq = BaseUtilsMore.simpleReverseIupacCodeComplement(motifSeq);
		int start = negStrand? -(motif.length()-1) : 0;
		int end = negStrand? 0 : motif.length()-1;
		

	
		for(int i = start, index = 0; i <= end; i++, index++){
			GenomeLoc loc = ref.getGenomeLocParser().createGenomeLoc(ref.getLocus().getContig(), ref.getLocus().getStart()+i );
			if( !ref.getWindow().containsP(loc) )
				return false;
			
			ReferenceContext tmpRef = new ReferenceContext(ref.getGenomeLocParser(),loc, ref.getWindow(),ref.getBases());
			refBytes[index] = tmpRef.getBase();
			if( !BaseUtils.isRegularBase(refBytes[index]) || !BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(motifSeq[index], refBytes[index]))
				return false;
		}
		

		
		return true;
	}
	
	private boolean isInCgi(RefMetaDataTracker tracker){
	//	for (RODRecordList rods : tracker.getBoundRodTracks()) {

	//		for (GATKFeature vc_input : rods) {
	//			if (vc_input != null ) {
	//				if(vc_input.getUnderlyingObject() instanceof BEDFeature){
	//					return true;
	//				}
					
	//			}
	//		}
	//	}
		return !tracker.getValues(cgi).isEmpty();
	}
	
	private double firstEndFraction(AlignmentContext context){
		int firstEnd = 0; 
		for (PileupElement p :  context.getBasePileup()) {
			if(!p.getRead().getReadPairedFlag()){
				return 1.0;
			}
			else{
				if(p.getRead().getFirstOfPairFlag()){
					firstEnd++;
				}
			}
		}
		return (double)firstEnd/(double)context.size();
	}
	
	private int mismatchFraction(ReferenceContext ref, AlignmentContext context, boolean firstEnd){
		

		int mismatches = 0;
		for (PileupElement p : context.getBasePileup()) {
			if(!p.getRead().getReadPairedFlag()){
				if(!firstEnd){
					return 0;
				}
				else{
					mismatches += !BaseUtils.basesAreEqual(ref.getBase(), p.getBase()) && BaseUtilsMore.isBisulfiteMismatch(ref.getBase(),p.getBase(), p.getRead().getReadNegativeStrandFlag(),false) ? 1 : 0;
					//System.err.println(ref.getBase() + "\t" + p.getBase() + "\t" + p.getRead().getReadNegativeStrandFlag() + "\t" + BaseUtilsMore.isBisulfiteMismatch(ref.getBase(),p.getBase(), p.getRead().getReadNegativeStrandFlag(),false));
				}
					
			}
			else{
				if(p.getRead().getFirstOfPairFlag() && firstEnd){
					mismatches += !BaseUtils.basesAreEqual(ref.getBase(), p.getBase()) && BaseUtilsMore.isBisulfiteMismatch(ref.getBase(),p.getBase(), p.getRead().getReadNegativeStrandFlag(),false) ? 1 : 0;
				}
				else if(p.getRead().getSecondOfPairFlag() && !firstEnd){
					mismatches += !BaseUtils.basesAreEqual(ref.getBase(), p.getBase()) && BaseUtilsMore.isBisulfiteMismatch(ref.getBase(),p.getBase(), p.getRead().getReadNegativeStrandFlag(),true) ? 1 : 0;
				}
			}
		}
		return mismatches;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.TreeReducible#treeReduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Long treeReduce(Long lhs, Long rhs) {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Long reduce(Long value, Long sum) {
		// TODO Auto-generated method stub
		return null;
	}


}
