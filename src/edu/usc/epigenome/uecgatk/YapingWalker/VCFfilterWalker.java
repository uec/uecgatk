/**
 * 
 */
package edu.usc.epigenome.uecgatk.YapingWalker;

import java.io.File;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.VariantEvalUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.kohsuke.args4j.Option;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time May 6, 2012 10:24:45 PM
 * 
 */
//@By(DataSource.REFERENCE_ORDERED_DATA)
@By(DataSource.READS)
@Requires({DataSource.READS,DataSource.REFERENCE, DataSource.REFERENCE_ORDERED_DATA, DataSource.REFERENCE_BASES})
@PartitionBy(PartitionType.LOCUS)
public class VCFfilterWalker extends LocusWalker<Integer, Integer> implements
		TreeReducible<Integer> {

	/**
	 * 
	 */
	@Input(fullName="old_vcfs", shortName = "oldVcfs", doc="input vcf file", required=true)
	 public List<RodBinding<VariantContext>> oldVcfs;
	
	@Output(fullName="new_vcf", shortName = "newVcf", doc="filtered vcf file", required=true)
	public VCFWriter newVcf = null;
	
	
	@Argument(shortName="minCov",doc="minimum covergae required for the position in VCF file, default:1", required=false)
	protected int minCov = 1;
	@Argument(shortName="maxCov",doc="maximum covergae required for the position in VCF file, default:120", required=false)
	protected int maxCov = 120;
	@Argument(shortName="minMapQ",doc="Minimum mapping quality in the covered loci when counting coverage, default:30", required=false)
	protected int minMapQ = 30;
	@Argument(shortName="minBaseQ",doc="Minimum base quality in the covered loci when counting coverage, default:0", required=false)
	protected int minBaseQ = 0;
	@Argument(shortName="minPercentOneStrand",doc="Minimum percentage of one strand only when total sequence coverage > 10X, default:0.2", required=false)
	protected double minPercentOneStrand = 0.2;
	@Argument(shortName="minNumOneStrand",doc="Minimum number of one strand only when less than 10 coverages, default:2", required=false)
	protected int minNumOneStrand = 2;
	@Argument(shortName="filterSB",doc="filter out loci that have high Strand bias", required=false)
	protected boolean filterSB = false;



//	private int recCounter = 1;
//	private int oldRecord = 0;
//	private int newRecord = 0;
	
	public void initialize(){
		Map<String, VCFHeader> headers = VCFUtils.getVCFHeadersFromRods(getToolkit(), oldVcfs);
		for(VCFHeader header : headers.values()){
			newVcf.writeHeader(header);
		}
	}
	
	
	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.LocusWalker#map(org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker, org.broadinstitute.sting.gatk.contexts.ReferenceContext, org.broadinstitute.sting.gatk.contexts.AlignmentContext)
	 */
	@Override
	public Integer map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		// TODO Auto-generated method stub
		if(tracker == null){
			return null;
		}
            // First, filter the VariantContext to represent only the samples for evaluation
		for(RodBinding<VariantContext> oldVcf : oldVcfs){
			for ( VariantContext vc_input : tracker.getValues(oldVcf, ref.getLocus()) ) {
				if ( vc_input != null ) { 
				//	oldRecord++;
					if(passFilter(context)){
	    			//	newRecord++;
						
						
	    				newVcf.add(vc_input);
	    			}
				}
			}
		}

		return null;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduceInit()
	 */
	@Override
	public Integer reduceInit() {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Integer reduce(Integer value, Integer sum) {
		// TODO Auto-generated method stub
		return null;
	}

	

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.TreeReducible#treeReduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Integer treeReduce(Integer lhs, Integer rhs) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public void onTraversalDone(Integer value) {
		//logger.info(String.format("Old VCF record: %d",oldRecord));
		//logger.info(String.format("New VCF record: %d",newRecord));
		logger.info("Finished!");
	}
	
	private boolean passFilter(AlignmentContext context){
		
		if(context.hasBasePileup()){
			int coverage = 0;
			int negStarndCov = 0;
			//System.err.println(context.hasReads() + "\t" + context.getPosition());
			for(PileupElement p : context.getBasePileup()){
				
				SAMRecord samRecord = p.getRead();
				
				if(samRecord.getDuplicateReadFlag() || samRecord.getNotPrimaryAlignmentFlag() || samRecord.getReadFailsVendorQualityCheckFlag() || samRecord.getReadUnmappedFlag() || samRecord.getMappingQuality() < minMapQ || p.getQual() < minBaseQ)
					continue;
				if(samRecord.getReadPairedFlag()){
					if (samRecord.getAlignmentStart() == samRecord.getMateAlignmentStart() && samRecord.getReadNegativeStrandFlag() == samRecord.getMateNegativeStrandFlag()){
						if (samRecord.getSecondOfPairFlag()) 
	 						continue;
					}
	 					if(!samRecord.getProperPairFlag())
	 						continue;
				}
				coverage++;
				if(samRecord.getReadNegativeStrandFlag()){
					negStarndCov++;
				}
			}
			if(coverage>=minCov && coverage<=maxCov){	
				if(filterSB){
					if((coverage>10 && (double)negStarndCov/(double)coverage > minPercentOneStrand && (double)(coverage-negStarndCov)/(double)coverage > minPercentOneStrand) || (coverage<=10 && negStarndCov > minNumOneStrand && (coverage-negStarndCov) > minNumOneStrand)){
						return true;
					}
				}
				else{
					return true;
				}
				
			}
		}
		return false;
	}
	
}
