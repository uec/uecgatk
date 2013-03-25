/**
 * 
 */
package edu.usc.epigenome.uecgatk.YapingWalker;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.broad.tribble.annotation.Strand;
import org.broad.tribble.bed.BEDFeature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
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

import edu.usc.epigenome.uecgatk.BisSNP.BaseUtilsMore;
import edu.usc.epigenome.uecgatk.BisSNP.BisSNPUtils;
import edu.usc.epigenome.uecgatk.YapingWalker.MotifFreqInGenomeWalker.Datum;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObject;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObjectWriterImp;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Feb 4, 2013 5:59:08 PM
 * 
 */
@Reference(window=@Window(start=-5000,stop=5000))
@Requires({DataSource.REFERENCE, DataSource.REFERENCE_BASES})
@By(DataSource.REFERENCE)
@Downsample(by=DownsampleType.NONE)
public class MotifPatternFeatureByBedWalker extends LocusWalker<MotifPatternFeatureByBedWalker.Datum, MotifPatternFeatureByBedWalker.Datum> implements TreeReducible<MotifPatternFeatureByBedWalker.Datum> {

	
	@Input(fullName = "motif", shortName = "motif", doc = "motif pattern to align", required = true)
	public String motif = null;
	
	@Input(fullName="segment", shortName = "seg", doc="segment bed file to align to", required=true) 
	 public RodBinding<BEDFeature> seg;
	
	@Output(fullName = "motif_file_name", shortName = "motifFile", doc = "Output file name with motif pattern frequency around the feature", required = true)
    public String motifFile = null;
	
	@Output(fullName = "seq_fwd_file_name", shortName = "seqFwdFile", doc = "Output file name with forward strand sequence composition around the feature", required = false)
    public String seqFwdFile = null;
	
	@Output(fullName = "seq_rev_file_name", shortName = "seqRevFile", doc = "Output file name with reverse strand sequence composition around the feature", required = false)
    public String seqRevFile = null;
	
	@Argument(fullName = "search_distance_to_feature", shortName = "distance", doc = "define the distance before or after feature, default: 200", required = false)
    public int distance = 200;
	
	@Argument(fullName = "bin_size", shortName = "binSize", doc = "define the bin size when sliding window. default: 1", required = false)
    public int binSize = 1;
	
	@Argument(fullName = "enable_orientation", shortName = "orientated", doc = "bed file are orientated by strand or not, default: not orientated", required = false)
    public boolean orientated = false;
	
	@Argument(fullName = "motif_alignment_type", shortName = "alignmentType", doc = "motif aligned at FiveEnd, ThreeEnd or Center of the feature file[FiveEnd, ThreeEnd, Center], default: Center", required = false)
    public MotifAlignmentType alignmentType = MotifAlignmentType.Center;
	
	@Argument(fullName = "ref_strand_to_align", shortName = "strand", doc = "which strand of reference genome are used in alignment[NONE, POSITIVE, NEGATIVE, BOTH], NONE means use bed defined strand, default: POSITIVE. ther is bug for both and negative..", required = false)
    public ReferenceStrandToUse strand = ReferenceStrandToUse.POSITIVE;
	
	private bedObjectWriterImp motifWriter = null;
	private bedObjectWriterImp seqFwdWriter = null;
	private bedObjectWriterImp seqRevWriter = null;
	
	private LinkedList<Double> tmpDataPointList = null;
	private LinkedList<Character> tmpFwdSeqList = null;
	private LinkedList<Character> tmpRevSeqList = null;
	
	private String preChr = null;
	private int preStart = -1;
	private int preEnd = -1;

	
	public void initialize(){
		tmpDataPointList = new LinkedList<Double>();
		tmpFwdSeqList = new LinkedList<Character>();
		tmpRevSeqList = new LinkedList<Character>();
		
		File fn1 = null;
		fn1 = new File(motifFile);
		motifWriter = new bedObjectWriterImp(fn1);
		
		if(seqFwdFile != null){
			File fn2 = null;
			fn2 = new File(seqFwdFile);
			seqFwdWriter = new bedObjectWriterImp(fn2);
		}
		if(seqRevFile != null){
			File fn3 = null;
			fn3 = new File(seqRevFile);
			seqRevWriter = new bedObjectWriterImp(fn3);
		}
		
		
		 
	}
	
	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.LocusWalker#map(org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker, org.broadinstitute.sting.gatk.contexts.ReferenceContext, org.broadinstitute.sting.gatk.contexts.AlignmentContext)
	 */
	@Override
	public Datum map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
		Datum value = new Datum();
		value.refLoci = 1L;
		if (tracker == null)
			return null;
		List<BEDFeature> values = tracker.getValues(seg);
		if(!values.isEmpty()){
			if(preChr == null){
				preChr = values.get(0).getChr();
				preStart = values.get(0).getStart();
				preEnd = values.get(0).getEnd();
			}
			else{
				if(values.get(0).getChr().equalsIgnoreCase(preChr) && values.get(0).getStart() == preStart && values.get(0).getEnd() == preEnd){
					return null;
				}
				else{
					preChr = values.get(0).getChr();
					preStart = values.get(0).getStart();
					preEnd = values.get(0).getEnd();
				}
			}
			Integer alignCenter = null;
			Strand bedStrand = values.get(0).getStrand();
			if(!orientated){
				
				if(alignmentType == MotifAlignmentType.Center){
					alignCenter = (int)((values.get(0).getStart() + values.get(0).getEnd())/2);
					
				}
				else if(alignmentType == MotifAlignmentType.FiveEnd){
					alignCenter = values.get(0).getStart();
				}
				else if(alignmentType == MotifAlignmentType.ThreeEnd){
					alignCenter = values.get(0).getEnd();
				}
			}
			else{
				if(alignmentType == MotifAlignmentType.Center){
					alignCenter = (int)((values.get(0).getStart() + values.get(0).getEnd())/2);
					
				}
				else if(alignmentType == MotifAlignmentType.FiveEnd){
					if(values.get(0).getStrand() == Strand.NEGATIVE){
						alignCenter = values.get(0).getEnd();
					}
					else{
						alignCenter = values.get(0).getStart();
					}
					
				}
				else if(alignmentType == MotifAlignmentType.ThreeEnd){
					if(values.get(0).getStrand() == Strand.NEGATIVE){
						alignCenter = values.get(0).getStart();
					}
					else{
						
						alignCenter = values.get(0).getEnd();
					}
					
				}
			}
			addMotifToList(alignCenter, strand, bedStrand, tmpDataPointList, ref);
			bedObject bedLine = new bedObject(values.get(0).getChr(), values.get(0).getStart(), values.get(0).getEnd(), values.get(0).getStrand(), (List)tmpDataPointList); 
			System.out.println(values.get(0).getChr() + "\t"  + values.get(0).getStart() + "\t"  +values.get(0).getEnd() + "\t"  +values.get(0).getStrand() + "\t"  + tmpDataPointList.size());
			motifWriter.add(bedLine);
			tmpDataPointList.clear();
			if(seqFwdFile != null){
				bedObject bedLineFwd = new bedObject(values.get(0).getChr(), values.get(0).getStart(), values.get(0).getEnd(), values.get(0).getStrand(), (List)tmpFwdSeqList); 
				//System.err.println(values.get(0).getChr() + "\t"  + values.get(0).getStart() + "\t"  +values.get(0).getEnd() + "\t"  +values.get(0).getStrand() + "\t"  + tmpFwdSeqList.size());
				seqFwdWriter.add(bedLineFwd);
				tmpFwdSeqList.clear();
			}
			if(seqRevFile != null){
				bedObject bedLineRev = new bedObject(values.get(0).getChr(), values.get(0).getStart(), values.get(0).getEnd(), values.get(0).getStrand(), (List)tmpRevSeqList); 
				//System.err.println(values.get(0).getChr() + "\t"  + values.get(0).getStart() + "\t"  +values.get(0).getEnd() + "\t"  +values.get(0).getStrand() + "\t"  + tmpRevSeqList.size());
				seqRevWriter.add(bedLineRev);
				tmpRevSeqList.clear();
			}
		}
		
		
		//value.motifLoci = (checkPattern(motif, ref, false) || checkPattern(motif, ref, true)) ?1L:0L;
		 return value;
	}
	
	public void onTraversalDone(Datum result) {
	//	logger.info("Genome size: " + result.refLoci);
	//	logger.info("Motif Pattern " + motif + " number: " + result.motifStat);
	//	logger.info("Motif frequency (%) in genome: " + String.format("%.3f", 100*(double)result.motifStat/(double)result.refLoci));
	//	logger.info("Bed file size: " + result.lociInBed);
	//	logger.info("Motif Pattern " + motif + " number in provided bed file: " + result.motifStatInBed);
	//	logger.info("Motif frequency (%) in provided bed file: " + String.format("%.3f", 100*(double)result.motifStatInBed/(double)result.lociInBed));	
		motifWriter.close();
		if(seqFwdFile != null){
			seqFwdWriter.close();
		}
		if(seqRevFile != null){
			seqRevWriter.close();
		}
		
		logger.info("Finished!");
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduceInit()
	 */
	@Override
	public Datum reduceInit() {

		return new Datum();
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Datum reduce(Datum value, Datum result) {
		//result.motifLoci += value.motifLoci;
	//	result.refLoci += value.refLoci;
	//	result.lociInBed += value.lociInBed;
	//	result.motifStat += value.motifStat;
	//	result.motifStatInBed += value.motifStatInBed;

		return result;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.TreeReducible#treeReduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Datum treeReduce(Datum lhs, Datum rhs) {
	//	lhs.refLoci += rhs.refLoci;
	//	lhs.lociInBed += rhs.lociInBed;
	//	lhs.motifStat += rhs.motifStat;
	//	lhs.motifStatInBed += rhs.motifStatInBed;
		
		return lhs;
	}
	
	
	
	public class Datum{
		public Long refLoci=0L;
		public Long lociInBed = 0L;
		//public Long motifLoci=0L;
		public Long motifStat = 0L;
		public Long motifStatInBed = 0L;
		

	}
	
	private void addMotifToList(int alignCenter, ReferenceStrandToUse refStrand, Strand bedStrand, LinkedList<Double> list, ReferenceContext ref){
		int start=alignCenter-distance;
		int end=alignCenter+distance;
		//if(orientated && bedStrand == Strand.NEGATIVE){
		//	start=alignCenter+distance;
		//	end=alignCenter-distance;
		//}
		for(int i=start; i<=end; i++){
			GenomeLoc loc = ref.getGenomeLocParser().createGenomeLoc(ref.getLocus().getContig(), i);
			if (!ref.getWindow().containsP(loc)) {
				if(orientated && bedStrand == Strand.NEGATIVE){
					list.offerFirst(Double.NaN);
					if(seqFwdFile != null){
						tmpFwdSeqList.offerFirst('N');
					}
					if(seqRevFile != null){
						tmpRevSeqList.offerFirst('N');
					}
				}
				else{
					list.offerLast(Double.NaN);
					if(seqFwdFile != null){
						tmpFwdSeqList.offerLast('N');
					}
					if(seqRevFile != null){
						tmpRevSeqList.offerLast('N');
					}
				}
				continue;
			}
			ReferenceContext tmpRef = new ReferenceContext(ref.getGenomeLocParser(), loc, ref.getWindow(), ref.getBases());
			byte base = tmpRef.getBase();
			int isPattern = -1;
			if(refStrand == ReferenceStrandToUse.NONE){
				if(bedStrand == Strand.NEGATIVE){
					isPattern = checkPattern(motif, tmpRef, true);
				}
				else{
					isPattern = checkPattern(motif, tmpRef, false);
				}
			}
			else if(refStrand == ReferenceStrandToUse.BOTH){
				isPattern = checkPattern(motif, tmpRef, false);
				if(isPattern == 0){
					isPattern=checkPattern(motif, tmpRef, true);
				} 
			}
			else if(refStrand == ReferenceStrandToUse.POSITIVE){
				isPattern = checkPattern(motif, tmpRef, false);
			}
			else if(refStrand == ReferenceStrandToUse.NEGATIVE){
				isPattern = checkPattern(motif, tmpRef, true);
			}
			
			if(orientated && bedStrand == Strand.NEGATIVE){
				list.offerFirst(isPattern==-1 ? Double.NaN : isPattern);
				if(seqFwdFile != null){
					tmpFwdSeqList.offerFirst((char)base);
				}
				if(seqRevFile != null){
					tmpRevSeqList.offerFirst((char)BaseUtils.simpleComplement(base));
				}
			}
			else{
				list.offerLast(isPattern==-1 ? Double.NaN : isPattern);
				if(seqFwdFile != null){
					tmpFwdSeqList.offerLast((char)base);
				}
				if(seqRevFile != null){
					tmpRevSeqList.offerLast((char)BaseUtils.simpleComplement(base));
				}
			}
			
		}
	}
	
	private int checkPattern(String motif, ReferenceContext ref, boolean negStrand){
		byte[] refBytes = new byte[motif.length()];
		byte[] motifSeq = motif.getBytes();
		if(negStrand)
			motifSeq = BaseUtilsMore.simpleIupacCodeComplement(motifSeq);
		int start = negStrand? -(motif.length()-1) : 0;
		int end = negStrand? 0 : motif.length()-1;
		
	
		for(int i = start, index = 0; i <= end; i++, index++){
			GenomeLoc loc = ref.getGenomeLocParser().createGenomeLoc(ref.getLocus().getContig(), ref.getLocus().getStart()+i );
			if( !ref.getWindow().containsP(loc) )
				return -1;
			
			ReferenceContext tmpRef = new ReferenceContext(ref.getGenomeLocParser(),loc, ref.getWindow(),ref.getBases());
			refBytes[index] = tmpRef.getBase();
			if( !BaseUtils.isRegularBase(refBytes[index]))
				return -1;
			if(!BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(motifSeq[index], refBytes[index]))
				return 0;
		}
		

		
		return 1;
	}

	private enum MotifAlignmentType{
		FiveEnd,
		ThreeEnd,
		Center
	}
	
	private enum ReferenceStrandToUse{
		POSITIVE,
		NEGATIVE,
		BOTH,
		NONE
	}

	
}
