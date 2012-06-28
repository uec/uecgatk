/**
 * maybe it is better to write rodWlker rather than Locus walker for this kind of alignment..
 */
package edu.usc.epigenome.uecgatk.NOMeSeqWalker;

import java.io.File;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.broad.tribble.Feature;
import org.broad.tribble.annotation.Strand;
import org.broad.tribble.bed.BEDFeature;
import org.broad.tribble.bed.FullBEDFeature;
import org.broad.tribble.bed.SimpleBEDFeature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.BAQMode;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Downsample;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteArgumentCollection;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteEnums;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteGenotyperEngine;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteVariantCallContext;
import edu.usc.epigenome.uecgatk.YapingWalker.NDRCallContext;
import edu.usc.epigenome.uecgatk.YapingWalker.NDRdetectWalker.windowsObject;
import edu.usc.epigenome.uecgatk.YapingWriter.GetCytosineContext.contextStatus;
import edu.usc.epigenome.uecgatk.YapingWriter.SortingBedObjectWriter;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObject;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObjectWriterImp;

import edu.usc.epigenome.uecgatk.YapingWriter.GetCytosineContext;
import edu.usc.epigenome.uecgatk.bisulfiteIndels.BisBAQ;
import edu.usc.epigenome.uecgatk.bisulfiteIndels.BisBAQMode;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 */

// require input bam file, reference sequence, dbSNP rod, -L interval list bed file, -B feature list bed file
@BisBAQMode(QualityMode = BisBAQ.QualityMode.ADD_TAG, ApplicationTime = BisBAQ.ApplicationTime.ON_INPUT)
@ReadFilters( {UnmappedReadFilter.class, BadMateFilter.class, NotPrimaryAlignmentFilter.class, DuplicateReadFilter.class} ) // Filter out all reads with zero mapping quality
@Reference(window=@Window(start=-500,stop=500))
@By(DataSource.REFERENCE)
@Downsample(by=DownsampleType.NONE)
public class MethyPatternFeatureWalker extends LocusWalker<Boolean, Boolean>
		implements TreeReducible<Boolean> {

	@ArgumentCollection private static BisulfiteArgumentCollection BAC = new BisulfiteArgumentCollection();
	
	@Input(fullName="aligned_feature", shortName = "feature" , doc="Input feature location to align", required=false)
	public RodBinding<BEDFeature> feature;
	
	@Output(fullName = "gch_file_name", shortName = "gchFile", doc = "Output GCH files name", required = true)
    public String gchFile = null;
	
	@Output(fullName = "wcg_file_name", shortName = "wcgFile", doc = "Output WCG files name", required = true)
    public String wcgFile = null;
	
	@Output(fullName = "hcg_file_name", shortName = "hcgFile", doc = "Output HCG files name", required = true)
    public String hcgFile = null;
	
	@Output(fullName = "gch_datapoint_file_name", shortName = "gchDataPoint", doc = "Output GCH data point files name", required = false)
    public String gchDataPoint = null;
	
//	@Argument(fullName = "feature_name", shortName = "feature", doc = "Feature name provide in -B:<name>,<type> <filename> option", required = false)
 //   public String feature = null;
	
	@Argument(fullName = "search_distance_to_feature", shortName = "distance", doc = "define the distance before or after feature", required = false)
    public int distance = 2000;
	
	@Argument(fullName = "bin_size", shortName = "binSize", doc = "define the bin size when sliding window. default: 1", required = false)
    public int binSize = 1;
	
	@Argument(fullName = "minium_CT_reads_count", shortName = "minCTdepth", doc = "minium number of CT reads should contained to calculate methylation value", required = false)
    public int minCTdepth = 1;
	
	//@Argument(fullName = "space_before_feature", shortName = "before", doc = "define the space before feature to detect", required = false)
  //  public int before = 2000;
	
	///@Argument(fullName = "space_after_feature", shortName = "after", doc = "define the space after feature to detect", required = false)
   // public int after = 2000;
	
	@Argument(fullName = "enable_orientation", shortName = "orientated", doc = "orientated by strand or not", required = false)
    public boolean orientated = false;
	
	@Argument(fullName = "motif_alignment_type", shortName = "alignmentType", doc = "motif aligned at FiveEnd, ThreeEnd or Center", required = false)
    public MotifAlignmentType alignmentType = MotifAlignmentType.Center;
	
	//output three files matrix file, contained GCH, WCG and HCG
	private bedObjectWriterImp gchWriter = null;
	private bedObjectWriterImp wcgWriter = null;
	private bedObjectWriterImp hcgWriter = null;
	
	private bedObjectWriterImp gchDataWriter = null;
	
	private boolean inFeature = false;
	private boolean writtenObject = false;
	
	private LinkedList<Double> tmpMethyValueListGch = null;
	private LinkedList<Double> tmpMethyValueListWcg = null;
	private LinkedList<Double> tmpMethyValueListHcg = null;
	
	private LinkedList<Integer> tmpDataPointListGch = null;
	
	private ReferenceOrderedDataSource rodIt = null;
	
	private BisulfiteGenotyperEngine BG_engine = null;
	
	private String chr = null;
	private int bedStart = 0;
	private int bedEnd = 0;
	private SimpleBEDFeature bed = null;
	private Strand strand = Strand.NONE;
	private FeatureCondition summary= null;
	
	public void initialize(){
		 File fn1 = new File(gchFile);
		 File fn2 = new File(wcgFile);
		 File fn3 = new File(hcgFile);
		 File fn4= null;
		 if(gchDataPoint != null)
			fn4 = new File(gchDataPoint);
		 gchWriter = new bedObjectWriterImp(fn1);
		 wcgWriter = new bedObjectWriterImp(fn2);
		 hcgWriter = new bedObjectWriterImp(fn3);
		 if(gchDataPoint != null)
			 gchDataWriter = new bedObjectWriterImp(fn4);
		 
		 tmpMethyValueListGch = new LinkedList<Double>();
		 tmpMethyValueListWcg = new LinkedList<Double>();
		 tmpMethyValueListHcg = new LinkedList<Double>();
		 tmpDataPointListGch = new LinkedList<Integer>();
		 rodIt = getToolkit().getRodDataSources().get(1);
		 BAC.sequencingMode = BisulfiteEnums.MethylSNPModel.GM;
		 BAC.makeCytosine();
		 summary = new FeatureCondition();
	}
	   
	@Override
	public Boolean map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {

	     GenomeLoc loc = ref.getLocus().getLocation();
    	 GenomeLoc searchLoc = getToolkit().getGenomeLocParser().createGenomeLoc(ref.getLocus().getLocation().getContig(), ref.getLocus().getLocation().getStart()-distance-10, ref.getLocus().getLocation().getStart()+distance+10);
    	 
    	 LocationAwareSeekableRODIterator locRodIt = rodIt.seek(searchLoc);
    	 if(!locRodIt.hasNext()){
    		 inFeature = false; 
    		// locRodIt.close();
    		 rodIt.close(locRodIt);
    		
    		 return null;
    	 }
    	 else{
    		   
    		 if(!inFeature){
    			 RODRecordList rodList = locRodIt.seekForward(searchLoc);
    			// System.err.println(rodList.get(0).getUnderlyingObject());
    			 //rodList.get(0).getUnderlyingObject()
    			// for(Object it : rodList){
    				 if(rodList.get(0).getUnderlyingObject() instanceof SimpleBEDFeature){
    					 SimpleBEDFeature bedTmp = (SimpleBEDFeature)rodList.get(0).getUnderlyingObject();
    					// System.err.println(rodList.get(0).getUnderlyingObject());
        				 if(loc.distance(getToolkit().getGenomeLocParser().createGenomeLoc(bedTmp.getChr(), (bedTmp.getStart() + bedTmp.getEnd())/2, (bedTmp.getStart() + bedTmp.getEnd())/2)) <= distance){
                			 bed = bedTmp;
                			 strand = bedTmp.getStrand();
            	    		 chr = bedTmp.getChr();
            	    		 bedStart = bedTmp.getStart();
            	    		 bedEnd = bedTmp.getEnd();
            	    		 inFeature = true;
            	    		 writtenObject = false;
            	    		// System.err.println(bed.getStart());
            	    		// break;
            	    	 }
        			// }
        			 
        			 }
        			 
    			 //} 
    		}
    		else{
     			// System.err.println(loc.distance(getToolkit().getGenomeLocParser().createGenomeLoc(bed.getChr(), (bed.getStart() + bed.getEnd())/2, (bed.getStart() + bed.getEnd())/2)) + "\t" + bed.toString());
     				 if(bed != null && loc.distance(getToolkit().getGenomeLocParser().createGenomeLoc(bed.getChr(), (bed.getStart() + bed.getEnd())/2, (bed.getStart() + bed.getEnd())/2)) > distance){

     	    		 inFeature = false;
     	    		 writtenObject = false;
     	    		// System.err.println(bed.getStart());
     	    		// break;
     				 }
     		}
    			 
            
    	 }
    	 
    	
    	 rodIt.close(locRodIt);
    	
	     if(!inFeature){
	    	 if(!writtenObject && (!tmpMethyValueListGch.isEmpty() || !tmpMethyValueListWcg.isEmpty() || !tmpMethyValueListHcg.isEmpty())){
	    		// System.err.println(tmpMethyValueListGch.size() + "\t" + tmpMethyValueListWcg.size());
	    		 if(tmpMethyValueListGch.size() == distance * 2 + 1){
	    			 
	    			 bedObject bedLineGch = new bedObject(chr, bedStart, bedEnd, strand, (List)tmpMethyValueListGch); //chr, start, end, strand, aveMethyNDR, gchNumNDR, gchDepthNDR, gchCTdepthNDR, aveMethyLinker, gchNumLinker, gchDepthLinker, gchCTdepthLinker
		    		 gchWriter.add(bedLineGch);
		    		 if(gchDataPoint != null){
		    			 bedObject bedLineGch2 = new bedObject(chr, bedStart, bedEnd, strand, (List)tmpDataPointListGch);
			    		 gchDataWriter.add(bedLineGch2);
		    		 }
		    		 
	    		 }
	    		
	    		 if(tmpMethyValueListWcg.size() == distance * 2 + 1){
	    			 bedObject bedLineWcg = new bedObject(chr, bedStart, bedEnd, strand, (List)tmpMethyValueListWcg); //chr, start, end, aveMethyNDR, gchNumNDR, gchDepthNDR, gchCTdepthNDR, aveMethyLinker, gchNumLinker, gchDepthLinker, gchCTdepthLinker
		    		 wcgWriter.add(bedLineWcg);
	    		 }
	    		 
	    		 if(tmpMethyValueListHcg.size() == distance * 2 + 1){
	    			 bedObject bedLineHcg = new bedObject(chr, bedStart, bedEnd, strand, (List)tmpMethyValueListHcg); //chr, start, end, aveMethyNDR, gchNumNDR, gchDepthNDR, gchCTdepthNDR, aveMethyLinker, gchNumLinker, gchDepthLinker, gchCTdepthLinker
		    		 hcgWriter.add(bedLineHcg);
	    		 }

	    		 tmpMethyValueListGch.clear();
	    		 tmpMethyValueListWcg.clear();
	    		 tmpMethyValueListHcg.clear();
	    		 if(gchDataPoint != null)
	    			 tmpDataPointListGch.clear();
	    		 
	    		 writtenObject = true;
	    		 logger.info(chr + "\t" + bedStart + "\t" + bedEnd);
	    	 }
	    	
	    	 return null;
	     }
	     else{

	 		if(context == null){
	 			
	 				addContextToList(null, strand, tmpMethyValueListGch);
	 				addContextToList(null, strand, tmpMethyValueListWcg);
	 				addContextToList(null, strand, tmpMethyValueListHcg);
	 				if(gchDataPoint != null)
	 					addCoverageToList(null, strand, tmpDataPointListGch);
	 			
	 				
	 			return null;
	 		}
	 		boolean isGch =false;
	 		boolean isWcg = false;
	 		boolean isHcg = false;
	 		BG_engine = new BisulfiteGenotyperEngine(tracker, ref, context, BAC, getToolkit());
	 		BisulfiteVariantCallContext bvc = BG_engine.getBisulfiteVariantCallContext();
	 		if(bvc == null || bvc.getSummaryAcrossRG().cytosinePatternConfirmedSet ==null){
	 			addContextToList(null, strand, tmpMethyValueListGch);
 				addContextToList(null, strand, tmpMethyValueListWcg);
 				addContextToList(null, strand, tmpMethyValueListHcg);
 				if(gchDataPoint != null)
 					addCoverageToList(null, strand, tmpDataPointListGch);
 				return null;
	 			
	 		}
	 		//	return null;
	 		HashSet<String> cytosinePatternConfirmedList = bvc.getSummaryAcrossRG().cytosinePatternConfirmedSet;
	 		
	 		for(String cytosinePattern : cytosinePatternConfirmedList){
	 			if(cytosinePattern.equalsIgnoreCase("GCH")){
	 				isGch=true;
	 			}
	 			else{
	 				if(cytosinePattern.equalsIgnoreCase("HCG")){
	 					isHcg=true;
		 			}
		 			if(cytosinePattern.equalsIgnoreCase("WCG")){
		 				isWcg=true;
		 			}
	 			}
	 				
	 		}
	 				
	 		if(isGch){
	 			
	 			addContextToList(bvc, strand, tmpMethyValueListGch);
	 			if(gchDataPoint != null)
	 				addCoverageToList(bvc, strand, tmpDataPointListGch);
	 			//if positive strand, offerLast(), if negative strand, offerFirst()
	 			
	 		}
	 		else{
	 			addContextToList(null, strand, tmpMethyValueListGch);
	 			if(gchDataPoint != null)
	 				addCoverageToList(null, strand, tmpDataPointListGch);
	 		}
	 		if(isWcg){
	 			
	 			addContextToList(bvc, strand, tmpMethyValueListWcg);
	 		}
	 		else{
	 			addContextToList(null, strand, tmpMethyValueListWcg);
	 		}
	 		if(isHcg){
	 			
	 			addContextToList(bvc, strand, tmpMethyValueListHcg);
	 		}
	 		else{
	 			addContextToList(null, strand, tmpMethyValueListHcg);
	 		}
	 		
	 		
	 		return null;
	     }

		
	}

	@Override
	public Boolean reduceInit() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Boolean reduce(Boolean value, Boolean sum) {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public Boolean treeReduce(Boolean lhs, Boolean rhs) {
		return null;
	}
	/*
	public FeatureCondition treeReduce(Boolean lhs, Boolean rhs) {
		// TODO Auto-generated method stub
		lhs.visitedFeatures += rhs.visitedFeatures;
		lhs.nNoGchValueFeatures += rhs.nNoGchValueFeatures;
		lhs.nNoHcgValueFeatures += rhs.nNoHcgValueFeatures;
		lhs.nNoWcgValueFeatures += rhs.nNoWcgValueFeatures;
		
		lhs.nGchConfidantlyCalled += rhs.nGchConfidantlyCalled;
		lhs.nHcgConfidantlyCalled += rhs.nHcgConfidantlyCalled;
		lhs.nWcgConfidantlyCalled += rhs.nWcgConfidantlyCalled;
		
		lhs.sumCTReadsGchConfidantlyCalled += rhs.sumCTReadsGchConfidantlyCalled;
		lhs.sumCTReadsHcgConfidantlyCalled += rhs.sumCTReadsGchConfidantlyCalled;
		lhs.sumCTReadsWcgConfidantlyCalled += rhs.sumCTReadsGchConfidantlyCalled;
		
		lhs.sumReadsGchConfidantlyCalled += rhs.sumReadsGchConfidantlyCalled;
		lhs.sumReadsHcgConfidantlyCalled += rhs.sumReadsHcgConfidantlyCalled;
		lhs.sumReadsWcgConfidantlyCalled += rhs.sumReadsWcgConfidantlyCalled;
		
		lhs.sumMethyReadsGchConfidantlyCalled += rhs.sumMethyReadsGchConfidantlyCalled;
		lhs.sumMethyReadsHcgConfidantlyCalled += rhs.sumMethyReadsHcgConfidantlyCalled;
		lhs.sumMethyReadsWcgConfidantlyCalled += rhs.sumMethyReadsWcgConfidantlyCalled;
		
		return lhs;
	}
	*/
	public void onTraversalDone(Boolean result) {
		
		gchWriter.close();
		wcgWriter.close();
		hcgWriter.close();
		if(gchDataPoint != null)
			gchDataWriter.close();
		logger.info("Finished!");
	}
	
	private void addContextToList(BisulfiteVariantCallContext bvc, Strand strand, LinkedList<Double> list){
		if(orientated){
			if(strand == Strand.NEGATIVE){
				if(bvc == null){
					list.offerFirst(Double.NaN);
				}
				else{
					list.offerFirst((double)bvc.getSummaryAcrossRG().numC/(double)(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT));
				}
			}
			else{
				if(bvc == null){
					list.offerLast(Double.NaN);
				}
				else{
					list.offerLast((double)bvc.getSummaryAcrossRG().numC/(double)(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT));
				}
			}
			
		}
		else{
			if(bvc == null){
				list.offerLast(Double.NaN);
			}
			else{
				list.offerLast((double)bvc.getSummaryAcrossRG().numC/(double)(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT));
			}
		}
	}
	
	private void addCoverageToList(BisulfiteVariantCallContext bvc, Strand strand, LinkedList<Integer> list){
		if(orientated){
			if(strand == Strand.NEGATIVE){
				if(bvc == null){
					list.offerFirst(0);
				}
				else{
					list.offerFirst(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT);
				}
			}
			else{
				if(bvc == null){
					list.offerLast(0);
				}
				else{
					list.offerLast(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT);
				}
			}
			
		}
		else{
			if(bvc == null){
				list.offerLast(0);
			}
			else{
				list.offerLast(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT);
			}
		}
	}
	
	private boolean listStatistics(LinkedList<Double> list, String type){
		boolean filter= true;
		int numNaElement = 0;
		if(type.equalsIgnoreCase("GCH")){
			
		}
		else if(type.equalsIgnoreCase("HCG")){
			
		}
		else if(type.equalsIgnoreCase("WCG")){
			
		}
		
		return filter;
	}
	
	private enum MotifAlignmentType{
		FiveEnd,
		ThreeEnd,
		Center
	}

	private class FeatureCondition{
		public long visitedFeatures = 0;
		public long nNoHcgValueFeatures = 0; //number of feature with not even one value of HCG (all are NA..)
		public long nNoWcgValueFeatures = 0; //number of feature with not even one value of WCG
		public long nNoGchValueFeatures = 0; //number of feature with not even one value of GCH
		
		public long nHcgConfidantlyCalled = 0;
		public long nWcgConfidantlyCalled = 0;
		public long nGchConfidantlyCalled = 0;
		
		public long sumReadsHcgConfidantlyCalled = 0;
		public long sumReadsWcgConfidantlyCalled = 0;
		public long sumReadsGchConfidantlyCalled = 0;
		
		public long sumCTReadsHcgConfidantlyCalled = 0;
		public long sumCTReadsWcgConfidantlyCalled = 0;
		public long sumCTReadsGchConfidantlyCalled = 0;
		
		public double sumMethyReadsHcgConfidantlyCalled = 0;
		public double sumMethyReadsWcgConfidantlyCalled = 0;
		public double sumMethyReadsGchConfidantlyCalled = 0;
		
		public double averageOfReadsCoveredHcg(){return (double)sumReadsHcgConfidantlyCalled/(double)nHcgConfidantlyCalled;}
		public double averageOfReadsCoveredWcg(){return (double)sumReadsWcgConfidantlyCalled/(double)nWcgConfidantlyCalled;}
		public double averageOfReadsCoveredGch(){return (double)sumReadsGchConfidantlyCalled/(double)nGchConfidantlyCalled;}
		
		public double averageOfCTReadsCoveredHcg(){return (double)sumCTReadsHcgConfidantlyCalled/(double)nHcgConfidantlyCalled;}
		public double averageOfCTReadsCoveredWcg(){return (double)sumCTReadsWcgConfidantlyCalled/(double)nWcgConfidantlyCalled;}
		public double averageOfCTReadsCoveredGch(){return (double)sumCTReadsGchConfidantlyCalled/(double)nGchConfidantlyCalled;}
		
		public double averageOfMethyCoveredHcg(){return sumMethyReadsHcgConfidantlyCalled/(double)nHcgConfidantlyCalled;}
		public double averageOfMethyCoveredWcg(){return sumMethyReadsWcgConfidantlyCalled/(double)nWcgConfidantlyCalled;}
		public double averageOfMethyCoveredGch(){return sumMethyReadsGchConfidantlyCalled/(double)nGchConfidantlyCalled;}
		
	}
	

	
}
