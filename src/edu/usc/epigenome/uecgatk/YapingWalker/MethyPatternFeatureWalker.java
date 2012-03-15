/**
 * 
 */
package edu.usc.epigenome.uecgatk.YapingWalker;

import java.io.File;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.broad.tribble.annotation.Strand;
import org.broad.tribble.bed.FullBEDFeature;
import org.broad.tribble.bed.SimpleBEDFeature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityReadFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentReadFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
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

import edu.usc.epigenome.uecgatk.YapingWalker.NDRCallContext;
import edu.usc.epigenome.uecgatk.YapingWalker.NDRdetectWalker.windowsObject;
import edu.usc.epigenome.uecgatk.YapingWriter.GetCytosineContext.contextStatus;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObject;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObjectWriterImp;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisSNPUtils;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisulfiteArgumentCollection;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisulfiteDiploidSNPGenotypePriors;
import edu.usc.epigenome.uecgatk.YapingWriter.GetCytosineContext;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 */

// require input bam file, reference sequence, dbSNP rod, -L interval list bed file, -B feature list bed file
@BAQMode(QualityMode = BAQ.QualityMode.OVERWRITE_QUALS, ApplicationTime = BAQ.ApplicationTime.ON_INPUT)
@ReadFilters( {UnmappedReadFilter.class, BadMateFilter.class, NotPrimaryAlignmentReadFilter.class, DuplicateReadFilter.class} ) // Filter out all reads with zero mapping quality
@Reference(window=@Window(start=-500,stop=500))
@By(DataSource.REFERENCE)
@Downsample(by=DownsampleType.NONE)
public class MethyPatternFeatureWalker extends LocusWalker<Boolean, Boolean>
		implements TreeReducible<Boolean> {

	@ArgumentCollection private static BisulfiteArgumentCollection BAC = new BisulfiteArgumentCollection();
	
	@Argument(fullName = "gch_file_name", shortName = "gchFile", doc = "Output GCH files name", required = true)
    public String gchFile = null;
	
	@Argument(fullName = "wcg_file_name", shortName = "wcgFile", doc = "Output WCG files name", required = true)
    public String wcgFile = null;
	
	@Argument(fullName = "hcg_file_name", shortName = "hcgFile", doc = "Output HCG files name", required = true)
    public String hcgFile = null;
	
	@Argument(fullName = "feature_name", shortName = "feature", doc = "Feature name provide in -B:<name>,<type> <filename> option", required = false)
    public String feature = null;
	
	@Argument(fullName = "search_distance_to_feature", shortName = "distance", doc = "define the distance before or after feature", required = false)
    public int distance = 2000;
	
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
	
	private boolean inFeature = false;
	private boolean writtenObject = false;
	
	private LinkedList<Double> tmpMethyValueListGch = null;
	private LinkedList<Double> tmpMethyValueListWcg = null;
	private LinkedList<Double> tmpMethyValueListHcg = null;
	
	private ReferenceOrderedDataSource rodIt = null;
	
	private String chr = null;
	private int bedStart = 0;
	private int bedEnd = 0;
	private SimpleBEDFeature bed = null;
	private Strand strand = null;

	private GenotypePriors genotypePriors;
	
	public void initialize(){
		 File fn1 = new File(gchFile);
		 File fn2 = new File(wcgFile);
		 File fn3 = new File(hcgFile);
		 gchWriter = new bedObjectWriterImp(fn1);
		 wcgWriter = new bedObjectWriterImp(fn2);
		 hcgWriter = new bedObjectWriterImp(fn3);
		 genotypePriors = new BisulfiteDiploidSNPGenotypePriors();
		 tmpMethyValueListGch = new LinkedList<Double>();
		 tmpMethyValueListWcg = new LinkedList<Double>();
		 tmpMethyValueListHcg = new LinkedList<Double>();
		 rodIt = getToolkit().getRodDataSources().get(0);
	}
	
	@Override
	public Boolean map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		// TODO Auto-generated method stub
		// List<Object> featureList = tracker.getReferenceMetaData(feature);
		// System.err.println(tracker.hasROD(feature) + feature + "\t" + featureList.size());
	 //    if ( featureList.size() < 1 || ! (featureList.get(0) instanceof SimpleBEDFeature)) {
	    	
	  //  	 throw new UserException.MalformedFile(String.format("%s track isn't a properly formated CallableBases object!", feature));
	   //  }
	    // Iterator<Object> iter = featureList.iterator();
	     
	     
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

    			// if(rodList.size() > 1){
    			//	 
    			// }
    		//	 else{
    				 SimpleBEDFeature bedTmp = (SimpleBEDFeature)rodList.get(0).getUnderlyingObject();
            		 if(loc.distance(getToolkit().getGenomeLocParser().createGenomeLoc(bedTmp.getChr(), (bedTmp.getStart() + bedTmp.getEnd())/2, (bedTmp.getStart() + bedTmp.getEnd())/2)) <= distance){
            			 bed = bedTmp;
            			 strand = bedTmp.getStrand();
        	    		 chr = bedTmp.getChr();
        	    		 bedStart = bedTmp.getStart();
        	    		 bedEnd = bedTmp.getEnd();
        	    		 inFeature = true;
        	    		 writtenObject = false;
        	    	//	 System.err.println(bed.getStart());
        	    		// break;
        	    	 }
    			// }
    			 
    		 }
    		 else{
    			// System.err.println(loc.distance(getToolkit().getGenomeLocParser().createGenomeLoc(bed.getChr(), (bed.getStart() + bed.getEnd())/2, (bed.getStart() + bed.getEnd())/2)));
    			 if(loc.distance(getToolkit().getGenomeLocParser().createGenomeLoc(bed.getChr(), (bed.getStart() + bed.getEnd())/2, (bed.getStart() + bed.getEnd())/2)) > distance){

    	    		 inFeature = false;
    	    		 writtenObject = false;
    	    	//	 System.err.println(bed.getStart());
    	    		// break;
    	    	 }
    		 }
    	 }
    	 
    	// System.err.println(loc.toString());
    	// System.err.println(searchLoc.toString());
    	 //locRodIt.close();
    	 rodIt.close(locRodIt);
    	// System.err.println(getToolkit().getRodDataSources().get(0).seek(searchLoc).seekForward(searchLoc).getLocation().toString());
    	// System.err.println(getToolkit().getRodDataSources().get(0).seek(searchLoc).seekForward(searchLoc).get(0).getLocation().toString());
    	// System.err.println(getToolkit().getRodDataSources().get(0).seek(searchLoc).peekNextLocation().getLocation().toString());
    	 /*
	     while(iter.hasNext()){
	    	 SimpleBEDFeature bed = (SimpleBEDFeature)iter.next();
	    	 
	    	 int motifCenter = 0;
				if(alignmentType == MotifAlignmentType.Center){
					motifCenter = (bed.getEnd() + bed.getStart())/2;
				}
				else if(alignmentType == MotifAlignmentType.ThreeEnd){
					motifCenter = bed.getEnd();
				}
				else if(alignmentType ==MotifAlignmentType.FiveEnd){
					motifCenter = bed.getStart();
				}
				else{
					System.err.println("alignment type error!");
				}
	    
	    	 if(ref.getLocus().getLocation().distance(loc) <= distance) {
	    		 strand = bed.getStrand();
	    		 chr = bed.getChr();
	    		 bedStart = bed.getStart();
	    		 bedEnd = bed.getEnd();
	    		 inFeature = true;
	    		 writtenObject = false;
	    		 break;
	    	 }
	    	 else{
	    		 inFeature = false; 
	    	 }
	    	 
	     }
	     */
	     if(!inFeature){
	    	 if(!writtenObject && (!tmpMethyValueListGch.isEmpty() || !tmpMethyValueListWcg.isEmpty() || !tmpMethyValueListHcg.isEmpty())){
	    		 if(tmpMethyValueListGch.size() == distance * 2 + 1){
	    			 bedObject bedLineGch = new bedObject(chr, bedStart, bedEnd, (List)tmpMethyValueListGch); //chr, start, end, aveMethyNDR, gchNumNDR, gchDepthNDR, gchCTdepthNDR, aveMethyLinker, gchNumLinker, gchDepthLinker, gchCTdepthLinker
		    		 gchWriter.add(bedLineGch);
	    		 }
	    		
	    		 if(tmpMethyValueListWcg.size() == distance * 2 + 1){
	    			 bedObject bedLineWcg = new bedObject(chr, bedStart, bedEnd, (List)tmpMethyValueListWcg); //chr, start, end, aveMethyNDR, gchNumNDR, gchDepthNDR, gchCTdepthNDR, aveMethyLinker, gchNumLinker, gchDepthLinker, gchCTdepthLinker
		    		 wcgWriter.add(bedLineWcg);
	    		 }
	    		 
	    		 if(tmpMethyValueListHcg.size() == distance * 2 + 1){
	    			 bedObject bedLineHcg = new bedObject(chr, bedStart, bedEnd, (List)tmpMethyValueListHcg); //chr, start, end, aveMethyNDR, gchNumNDR, gchDepthNDR, gchCTdepthNDR, aveMethyLinker, gchNumLinker, gchDepthLinker, gchCTdepthLinker
		    		 hcgWriter.add(bedLineHcg);
	    		 }

	    		 tmpMethyValueListGch.clear();
	    		 tmpMethyValueListWcg.clear();
	    		 tmpMethyValueListHcg.clear();
	    		 writtenObject = true;
	    		 logger.info(chr + "\t" + bedStart + "\t" + bedEnd);
	    	 }
	    	
	    	 return null;
	     }
	     else{
	    	 String cytosinePatternGch = "GCH-2";
	 		String cytosinePatternWcg = "WCG-2";
	 		String cytosinePatternHcg = "HCG-2";
	 		double methyStatusGch = BAC.forceGch; //H1: 0.36; imr90:0.45
	 		double methyStatusCpg =BAC.forceCpg; //0.80

	 		
	 		BisSNPUtils it = new BisSNPUtils(BAC);
	 		AlignmentContext stratifiedContexts = it.getFilteredAndStratifiedContexts(BAC, ref, context);
	 		
	 		if(stratifiedContexts == null){
	 			
	 				addContextToList(null, strand, tmpMethyValueListGch);
	 				addContextToList(null, strand, tmpMethyValueListWcg);
	 				addContextToList(null, strand, tmpMethyValueListHcg);
	 			
	 				
	 			return null;
	 		}
	 		
	 		boolean isGch = it.checkCytosineStatus(cytosinePatternGch, stratifiedContexts.getBasePileup(), tracker, ref, (BisulfiteDiploidSNPGenotypePriors) genotypePriors, BAC, methyStatusGch);
	 		boolean isWcg = it.checkCytosineStatus(cytosinePatternWcg, stratifiedContexts.getBasePileup(), tracker, ref, (BisulfiteDiploidSNPGenotypePriors) genotypePriors, BAC, methyStatusCpg);
	 		boolean isHcg = it.checkCytosineStatus(cytosinePatternHcg, stratifiedContexts.getBasePileup(), tracker, ref, (BisulfiteDiploidSNPGenotypePriors) genotypePriors, BAC, methyStatusCpg);
	 		GetCytosineContext contextToTest = new GetCytosineContext();
	 		boolean paired = stratifiedContexts.getBasePileup().getReads().get(0).getReadPairedFlag();
 			contextStatus status = contextToTest.getContext(stratifiedContexts.getBasePileup(), paired, minCTdepth, minCTdepth);		
	 		if(isGch){
	 			
	 			addContextToList(status, strand, tmpMethyValueListGch);
	 			//if positive strand, offerLast(), if negative strand, offerFirst()
	 			
	 		}
	 		else{
	 			addContextToList(null, strand, tmpMethyValueListGch);
	 		}
	 		if(isWcg){
	 			
	 			addContextToList(status, strand, tmpMethyValueListWcg);
	 		}
	 		else{
	 			addContextToList(null, strand, tmpMethyValueListWcg);
	 		}
	 		if(isHcg){
	 			
	 			addContextToList(status, strand, tmpMethyValueListHcg);
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
		// TODO Auto-generated method stub
		return null;
	}
	
	public void onTraversalDone(Boolean over) {
		
		gchWriter.close();
		wcgWriter.close();
		hcgWriter.close();
		logger.info("Finished!");
	}
	
	private void addContextToList(contextStatus status, Strand strand, LinkedList<Double> list){
		if(orientated){
			if(strand == Strand.NEGATIVE){
				if(status == null){
					list.offerFirst(Double.NaN);
				}
				else{
					list.offerFirst((double)status.numC/(double)(status.numC + status.numT));
				}
			}
			else{
				if(status == null){
					list.offerLast(Double.NaN);
				}
				else{
					list.offerLast((double)status.numC/(double)(status.numC + status.numT));
				}
			}
			
		}
		else{
			if(status == null){
				list.offerLast(Double.NaN);
			}
			else{
				list.offerLast((double)status.numC/(double)(status.numC + status.numT));
			}
		}
	}
	
	private enum MotifAlignmentType{
		FiveEnd,
		ThreeEnd,
		Center
	}

}
