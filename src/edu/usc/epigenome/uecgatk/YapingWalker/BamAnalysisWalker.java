/**
 * 
 */
package edu.usc.epigenome.uecgatk.YapingWalker;

import java.io.File;
import java.io.PrintStream;
import java.util.LinkedList;
import java.util.List;

import org.broad.tribble.bed.BEDFeature;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Downsample;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import edu.usc.epigenome.uecgatk.YapingWriter.SortingBedObjectWriter;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObject;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObjectWriterImp;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time May 18, 2012 12:24:10 PM
 * 
 */
@ReadFilters( {} ) // Filter out all reads with zero mapping quality
@Reference(window=@Window(start=-200,stop=200))
@By(DataSource.REFERENCE)
@Downsample(by=DownsampleType.NONE)
public class BamAnalysisWalker extends LocusWalker<Long, Long> implements
		TreeReducible<Long> {

	/**
	 * 
	 */
	@Input(fullName="cpg_island_feature", shortName = "cgi" , doc="Input cpg island feature location", required=false)
	public RodBinding<BEDFeature> cgi;
	
	@Output(doc = "Output summary files", required = true)
    public String file = null;
	
	private bedObjectWriterImp writer = null;
	private SortingBedObjectWriter multiThreadWriter = null;
	private static int MAXIMUM_CACHE_FOR_OUTPUT = 100000;
	
	private double IN_ELEMENT = 1.0;
	private double NO_ELEMENT = 0.0;
	private double NO_VALUE = Double.NaN;
	
	public void initialize(){
		writer = new bedObjectWriterImp(new File(file));
		if(getToolkit().getArguments().numberOfThreads > 1){
			multiThreadWriter = new SortingBedObjectWriter(writer, MAXIMUM_CACHE_FOR_OUTPUT);
		}
		String header = "chr\tstart\tend\tstrand\tcoverage\tcgi\tref_cg\tc_reads\tt_reads\t1st_end_reads\t2nd_end_reads\tfwd_strand\trev_strand\tproper_paired\t" +
				"not_primary_aligned\tunmapped\tduplicated\tfailed_vender_checked\tinverse_duplicated\tpass_filter_coverage\tpass_filter_c_reads\tpass_filter_t_reads\tpass_filter_1st_end_reads\tpass_filter_2nd_end_reads\tpass_filter_fwd_reads\tpass_filter_rev_reads\n";
		if(getToolkit().getArguments().numberOfThreads > 1){
			multiThreadWriter.writeHeader(header);
		}
		else{
			writer.addHeader(header);
		}
		
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.LocusWalker#map(org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker, org.broadinstitute.sting.gatk.contexts.ReferenceContext, org.broadinstitute.sting.gatk.contexts.AlignmentContext)
	 */
	@Override //todo: should seperate it into different read group
	public Long map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		LinkedList<Double> stats =  new LinkedList<Double>();
		String chr = ref.getLocus().getContig();
		int bedStart = ref.getLocus().getStart()-1;
		int bedEnd = ref.getLocus().getStart();
		stats.addLast((double)context.size());
		//if(ref.getLocus().getStart()==7000000)
		//System.err.println(context.size());
		if(tracker.hasValues(cgi) && tracker.getValues(cgi).get(0) instanceof BEDFeature){
			stats.addLast(IN_ELEMENT);
		}
		else{
			stats.addLast(NO_ELEMENT);
		}
		stats.addLast(BaseUtils.basesAreEqual(BaseUtils.C, ref.getBase()) || BaseUtils.basesAreEqual(BaseUtils.G, ref.getBase()) ?  IN_ELEMENT : NO_ELEMENT);	
		if(context.size() != 0){
			statsOnPileup(context.getBasePileup(),stats);
		}
		else{
			statsOnPileup(null,stats);
		}
		bedObject bedLine = new bedObject(chr, bedStart, bedEnd, (List)stats);
		if(getToolkit().getArguments().numberOfThreads > 1){
			multiThreadWriter.add(bedLine);
		}
		else{
			writer.add(bedLine);
		}
		
		return null;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduceInit()
	 */
	@Override
	public Long reduceInit() {
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
	
	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.TreeReducible#treeReduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Long treeReduce(Long lhs, Long rhs) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public void onTraversalDone(Boolean result) {
		if(getToolkit().getArguments().numberOfThreads > 1){
			multiThreadWriter.close();
		}
		else{
			writer.close();
		}
		writer.close();
		logger.info("Finished!");
	}

	private void statsOnPileup(ReadBackedPileup pileup, LinkedList<Double> list){
		if(pileup == null){
			fillList(list);
			return;
		}
		double cReads = 0;
		double tReads = 0;
		double firstEndReads = 0;
		double secondEndReads = 0;
		double fwdReads = 0;
		double revReads = 0;
		double properPairedReads = 0;
		double NotPrimaryAlignedReads = 0;
		double UnmappedReads = 0;
		double DuplicatedReads = 0;
		double FailedVenderCheckReads = 0;
		double InverseDuplicatedReads = 0;
		double PassFilterReads = 0;
		double PassFilterCReads = 0;
		double PassFilterTReads = 0;
		double PassFilterFirstEndReads = 0;
		double PassFilterSecondEndReads = 0;
		double PassFilterFwdReads = 0;
		double PassFilterRevReads = 0;
		
		for(PileupElement p : pileup){
			//if(pileup.getLocation().getStart()==7000000)
			//	System.err.println((char)p.getBase());
			if(p.getRead().getReadPairedFlag() && p.getRead().getSecondOfPairFlag()){
				if(p.getRead().getReadNegativeStrandFlag()){
					cReads += BaseUtils.basesAreEqual(BaseUtils.C, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
					tReads += BaseUtils.basesAreEqual(BaseUtils.T, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
				}
				else{
					cReads += BaseUtils.basesAreEqual(BaseUtils.G, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
					tReads += BaseUtils.basesAreEqual(BaseUtils.A, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
				}
				
			}
			else{
				if(p.getRead().getReadNegativeStrandFlag()){
					cReads += BaseUtils.basesAreEqual(BaseUtils.G, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
					tReads += BaseUtils.basesAreEqual(BaseUtils.A, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
				}
				else{
					cReads += BaseUtils.basesAreEqual(BaseUtils.C, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
					tReads += BaseUtils.basesAreEqual(BaseUtils.T, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
				}
			}
			
			firstEndReads += p.getRead().getReadPairedFlag() && p.getRead().getFirstOfPairFlag() ? IN_ELEMENT : NO_ELEMENT;
			secondEndReads += p.getRead().getReadPairedFlag() && p.getRead().getSecondOfPairFlag() ? IN_ELEMENT : NO_ELEMENT;
			fwdReads += p.getRead().getReadNegativeStrandFlag() ? NO_ELEMENT : IN_ELEMENT;
			revReads += p.getRead().getReadNegativeStrandFlag() ? IN_ELEMENT : NO_ELEMENT;
			
			properPairedReads += p.getRead().getReadPairedFlag() && p.getRead().getProperPairFlag() ? IN_ELEMENT : NO_ELEMENT;
			
			NotPrimaryAlignedReads += p.getRead().getNotPrimaryAlignmentFlag() ? IN_ELEMENT : NO_ELEMENT;
			UnmappedReads += p.getRead().getReadUnmappedFlag() ? IN_ELEMENT : NO_ELEMENT;
			DuplicatedReads += p.getRead().getDuplicateReadFlag() ? IN_ELEMENT : NO_ELEMENT;
			FailedVenderCheckReads += p.getRead().getReadFailsVendorQualityCheckFlag() ? IN_ELEMENT : NO_ELEMENT;
			InverseDuplicatedReads += p.getRead().getReadPairedFlag() && p.getRead().getAlignmentStart() == p.getRead().getMateAlignmentStart() && p.getRead().getReadNegativeStrandFlag() == p.getRead().getMateNegativeStrandFlag() ? IN_ELEMENT : NO_ELEMENT;
			if(!(p.getRead().getNotPrimaryAlignmentFlag() && p.getRead().getReadUnmappedFlag() && p.getRead().getDuplicateReadFlag() && p.getRead().getReadFailsVendorQualityCheckFlag())){
				if(p.getRead().getReadPairedFlag()){
					if(p.getRead().getProperPairFlag()){
						if(p.getRead().getFirstOfPairFlag()){
							PassFilterFirstEndReads++;
							PassFilterReads++;
							PassFilterFwdReads += p.getRead().getReadNegativeStrandFlag() ? NO_ELEMENT : IN_ELEMENT;
							PassFilterRevReads += p.getRead().getReadNegativeStrandFlag() ? IN_ELEMENT : NO_ELEMENT;
							if(p.getRead().getReadNegativeStrandFlag()){
								PassFilterCReads += BaseUtils.basesAreEqual(BaseUtils.G, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
								PassFilterTReads += BaseUtils.basesAreEqual(BaseUtils.A, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
							}
							else{
								PassFilterCReads += BaseUtils.basesAreEqual(BaseUtils.C, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
								PassFilterTReads += BaseUtils.basesAreEqual(BaseUtils.T, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
							}
							
						}
						else{
							if(!(p.getRead().getAlignmentStart() == p.getRead().getMateAlignmentStart() && p.getRead().getReadNegativeStrandFlag() == p.getRead().getMateNegativeStrandFlag())){
								PassFilterSecondEndReads++;
								PassFilterReads++;
								PassFilterFwdReads += p.getRead().getReadNegativeStrandFlag() ? NO_ELEMENT : IN_ELEMENT;
								PassFilterRevReads += p.getRead().getReadNegativeStrandFlag() ? IN_ELEMENT : NO_ELEMENT;
								if(p.getRead().getReadNegativeStrandFlag()){
									PassFilterCReads += BaseUtils.basesAreEqual(BaseUtils.C, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
									PassFilterTReads += BaseUtils.basesAreEqual(BaseUtils.T, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
								}
								else{
									PassFilterCReads += BaseUtils.basesAreEqual(BaseUtils.G, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
									PassFilterTReads += BaseUtils.basesAreEqual(BaseUtils.A, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
								}
							}
						}
					}
				}
				else{
					PassFilterFirstEndReads++;
					PassFilterReads++;
					PassFilterFwdReads += p.getRead().getReadNegativeStrandFlag() ? NO_ELEMENT : IN_ELEMENT;
					PassFilterRevReads += p.getRead().getReadNegativeStrandFlag() ? IN_ELEMENT : NO_ELEMENT;
					if(p.getRead().getReadNegativeStrandFlag()){
						PassFilterCReads += BaseUtils.basesAreEqual(BaseUtils.G, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
						PassFilterTReads += BaseUtils.basesAreEqual(BaseUtils.A, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
					}
					else{
						PassFilterCReads += BaseUtils.basesAreEqual(BaseUtils.C, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
						PassFilterTReads += BaseUtils.basesAreEqual(BaseUtils.T, p.getBase()) ? IN_ELEMENT : NO_ELEMENT;
					}
				}
			}
		}
		list.addLast(cReads);
		list.addLast(tReads);
		list.addLast(firstEndReads);
		list.addLast(secondEndReads);
		list.addLast(fwdReads);
		list.addLast(revReads);
		list.addLast(properPairedReads);
		list.addLast(NotPrimaryAlignedReads);
		list.addLast(UnmappedReads);
		list.addLast(DuplicatedReads);
		list.addLast(FailedVenderCheckReads);
		list.addLast(InverseDuplicatedReads);
		list.addLast(PassFilterReads);
		list.addLast(PassFilterCReads);
		list.addLast(PassFilterTReads);
		list.addLast(PassFilterFirstEndReads);
		list.addLast(PassFilterSecondEndReads);
		list.addLast(PassFilterFwdReads);
		list.addLast(PassFilterRevReads);
		
		
	}
	
	private void fillList(LinkedList<Double> list){
		for(int i = 1; i<=19; i++){
			list.addLast(NO_VALUE);
		}
	}
}
