/**
 * 
 */
package edu.usc.epigenome.uecgatk.YapingWalker;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMSequenceDictionary;

import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
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
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import edu.usc.epigenome.uecgatk.BisSNP.BaseUtilsMore;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jun 20, 2012 4:45:10 PM
 * 
 */


@Reference(window=@Window(start=-500,stop=500))
@Requires({DataSource.REFERENCE, DataSource.REFERENCE_BASES})
@By(DataSource.REFERENCE)
@Downsample(by=DownsampleType.NONE)
public class MotifFreqInGenomeWalker extends LocusWalker<MotifFreqInGenomeWalker.Datum, MotifFreqInGenomeWalker.Datum> implements TreeReducible<MotifFreqInGenomeWalker.Datum> {

	@Input(fullName = "motif_to_search", shortName = "motifs", doc = "motif pattern to search in provided genome", required = true)
	public ArrayList<String> motifs = null;
	
	private IndexedFastaSequenceFile referenceReader;
	/**
	 * 
	 */
	public void initialize(){
		try {
            referenceReader = new CachingIndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(getToolkit().getArguments().referenceFile,ex);
        }
	}

	

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.LocusWalker#map(org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker, org.broadinstitute.sting.gatk.contexts.ReferenceContext, org.broadinstitute.sting.gatk.contexts.AlignmentContext)
	 */
	@Override
	public Datum map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext contexts) {
		Datum value = new Datum(motifs);
		value.refLoci = 1L;
		for(String motif : motifs)
			value.motifStat.put(motif, (checkPattern(motif, ref, false) || checkPattern(motif, ref, true)) ?1L:0L);
		//value.motifLoci = (checkPattern(motif, ref, false) || checkPattern(motif, ref, true)) ?1L:0L;
		 return value;
	}
	

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduceInit()
	 */
	@Override
	public Datum reduceInit() {

		return new Datum(motifs);
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Datum reduce(Datum value, Datum result) {
		
		//result.motifLoci += value.motifLoci;
		result.refLoci += value.refLoci;
		for(String motif : motifs)
			result.motifStat.put(motif, result.motifStat.get(motif) + value.motifStat.get(motif));
		return result;
	}


	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.TreeReducible#treeReduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Datum treeReduce(Datum lhs, Datum rhs) {
		
		lhs.refLoci += rhs.refLoci;
		for(String motif : motifs)
			lhs.motifStat.put(motif, lhs.motifStat.get(motif) + rhs.motifStat.get(motif));
		return lhs;
	}
	
	public void onTraversalDone(Datum result) {
		logger.info("Genome size: " + result.refLoci);
		for(String motif : motifs){
			logger.info("Motif Pattern " + motif + " number: " + result.motifStat.get(motif));
			logger.info("Motif frequency (%): " + String.format("%.3f", 100*(double)result.motifStat.get(motif)/(double)result.refLoci));
		}
		
		logger.info("Finished!");
	}
	
	private boolean checkPattern(String motif, ReferenceContext ref, boolean negStrand){
		byte[] refBytes = new byte[motif.length()];
		byte[] motifSeq = motif.getBytes();
		if(negStrand)
			motifSeq = BaseUtilsMore.simpleReverseIupacCodeComplement(motifSeq);
		int start = negStrand? -(motif.length()-1) : 0;
		int end = negStrand? 0 : motif.length()-1;
		
	//	if(ref.getLocus().getStart()+start >= 0 && ref.getLocus().getStart()+end < referenceReader.getSequence(ref.getLocus().getContig()).length()){
	//		refBytes = referenceReader.getSubsequenceAt(ref.getLocus().getContig(), ref.getLocus().getStart()+start, ref.getLocus().getStart()+end).getBases();
	//	}
	//	else{
	//		return false;
	//	}
	
		for(int i = start, index = 0; i <= end; i++, index++){
			GenomeLoc loc = ref.getGenomeLocParser().createGenomeLoc(ref.getLocus().getContig(), ref.getLocus().getStart()+i );
			if( !ref.getWindow().containsP(loc) )
				return false;
			
			ReferenceContext tmpRef = new ReferenceContext(ref.getGenomeLocParser(),loc, ref.getWindow(),ref.getBases());
			refBytes[index] = tmpRef.getBase();
			if( !BaseUtils.isRegularBase(refBytes[index]) || !BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(motifSeq[index], refBytes[index]))
				return false;
		}
		
	//	System.err.println(ref.getLocus() + "\t" + negStrand);
	//	for(int i = start, index = 0; i <= end; i++, index++){
	//		if( !BaseUtils.isRegularBase(refBytes[index]) || !BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(motifSeq[index], refBytes[index]))
	//			return false;
	//	}
	//	System.err.println(new String(refBytes));
	//	System.err.println(new String(motifSeq));
	//	System.err.println(ref.getLocus() + "\t" + negStrand);
		
		return true;
	}
	
	public class Datum{
		public Long refLoci=0L;
		//public Long motifLoci=0L;
		public HashMap<String,Long> motifStat = null;
		public Datum(ArrayList<String> motifs){
			motifStat = new HashMap<String,Long>();
			for(String motif : motifs)
				motifStat.put(motif, 0L);
		}
	}
	
}
