/**
 * 
 */
package edu.usc.epigenome.uecgatk.YapingWalker;

import java.io.File;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.broad.tribble.annotation.Strand;
import org.broad.tribble.bed.FullBEDFeature;
import org.broad.tribble.bed.SimpleBEDFeature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.MappingQualityReadFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentReadFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
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
@ReadFilters( {UnmappedReadFilter.class, BadMateFilter.class, NotPrimaryAlignmentReadFilter.class, DuplicateReadFilter.class, MappingQualityReadFilter.class} ) // Filter out all reads with zero mapping quality
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
	
	@Argument(fullName = "feature_name", shortName = "feature", doc = "Feature name provide in -B <name>,<type>,<file> option", required = true)
    public String feature = null;
	
	@Argument(fullName = "search_distance_to_feature", shortName = "distance", doc = "define the distance before or after feature", required = false)
    public int distance = 2000;
	
	//@Argument(fullName = "space_before_feature", shortName = "before", doc = "define the space before feature to detect", required = false)
  //  public int before = 2000;
	
	///@Argument(fullName = "space_after_feature", shortName = "after", doc = "define the space after feature to detect", required = false)
   // public int after = 2000;
	
	@Argument(fullName = "enable_orientation", shortName = "orientated", doc = "orientated by strand or not", required = false)
    public boolean orientated = false;
	
	//output three files matrix file, contained GCH, WCG and HCG
	private bedObjectWriterImp gchWriter = null;
	private bedObjectWriterImp wcgWriter = null;
	private bedObjectWriterImp hcgWriter = null;
	
	private boolean inFeature = false;
	private boolean writtenObject = false;
	
	private LinkedList<Double> tmpMethyValueListGch = null;
	private LinkedList<Double> tmpMethyValueListWcg = null;
	private LinkedList<Double> tmpMethyValueListHcg = null;

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
	}
	
	@Override
	public Boolean map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		// TODO Auto-generated method stub
		 List<Object> featureList = tracker.getReferenceMetaData(feature);
	     if ( featureList.size() <= 1 || ! (featureList.get(0) instanceof SimpleBEDFeature)) {
	    	 throw new UserException.MalformedFile(String.format("%s track isn't a properly formated CallableBases object!", feature));
	     }
	     Iterator<Object> iter = featureList.iterator();
	     
	     while(iter.hasNext()){
	    	 SimpleBEDFeature bed = (SimpleBEDFeature)iter.next();
	    	 String contextChr = context.getContig();
	    	 int contextPos = context.getLocation().getStart();
	    	 String chr = context.getContig();
	    	 int bedPos = bed.getStart();
	    	 Strand strand = bed.getStrand();
	    	 if(chr.equalsIgnoreCase(contextChr)){
	    		 if(strand == Strand.POSITIVE){
		    		 
		    	 }
		    	 else{
		    		 
		    	 } 
	    	 }
	    	 
	    	 GenomeLoc loc = getToolkit().getGenomeLocParser().createGenomeLoc(bed.getChr(), bed.getStart(), bed.getEnd());
	    	 if(context.getLocation().distance(loc) <= distance) {
	    		 inFeature = true;
	    		 writtenObject = false;
	    		 break;
	    	 }
	    	 else{
	    		 inFeature = false; 
	    	 }
	    	 
	     }
	     if(!inFeature && writtenObject){
	    	 
	    	 return null;
	     }
	     
	     
	        
		String cytosinePatternGch = "GCH-2";
		String cytosinePatternWcg = "WCG-2";
		String cytosinePatternHcg = "HCG-2";
		double methyStatusGch = 0.36; //H1: 0.36; imr90:0.45
		double methyStatusCpg = 0.80; //H1: 0.36; imr90:0.45

		BisSNPUtils it = new BisSNPUtils(BAC);
		AlignmentContext stratifiedContexts = it.getFilteredAndStratifiedContexts(BAC, ref, context);
		
		if(it.checkCytosineStatus(cytosinePatternGch, stratifiedContexts.getBasePileup(), tracker, ref, (BisulfiteDiploidSNPGenotypePriors) genotypePriors, BAC, methyStatusGch)){
			GetCytosineContext contextToTest = new GetCytosineContext();
			contextStatus status = contextToTest.getContext(stratifiedContexts.getBasePileup(), BAC.pairedEndMode);
			//if positive strand, offerLast(), if negative strand, offerFirst()
			
		}
		if(it.checkCytosineStatus(cytosinePatternWcg, stratifiedContexts.getBasePileup(), tracker, ref, (BisulfiteDiploidSNPGenotypePriors) genotypePriors, BAC, methyStatusCpg)){
			
		}
		if(it.checkCytosineStatus(cytosinePatternHcg, stratifiedContexts.getBasePileup(), tracker, ref, (BisulfiteDiploidSNPGenotypePriors) genotypePriors, BAC, methyStatusCpg)){
			
		}
		
		
		
		return null;
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

}
