/**
 * 
 */
package edu.usc.epigenome.uecgatk.NOMeSeqWalker;

import java.io.File;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

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
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Downsample;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.GenomeLoc;

import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteArgumentCollection;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteGenotyperEngine;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteVariantCallContext;
import edu.usc.epigenome.uecgatk.YapingWalker.MotifFreqInGenomeWalker.Datum;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObject;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObjectWriterImp;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jun 25, 2012 4:52:28 PM
 * 
 */
@Reference(window=@Window(start=-500,stop=500))
@Requires({DataSource.REFERENCE, DataSource.REFERENCE_BASES})
@By(DataSource.REFERENCE)
@Downsample(by=DownsampleType.NONE)
public class MethyPatternFeatureByBedWalker extends LocusWalker<Boolean, Boolean> implements TreeReducible<Boolean> {

	
	@Input(fullName="input_values_for_align", shortName = "values" , doc="Input methylation values for alignment", required=true)
	public RodBinding<BEDFeature> values;

	@Input(fullName="aligned_feature", shortName = "feature" , doc="Input feature location to align", required=true)
	public RodBinding<BEDFeature> feature;
	
	@Output(fullName = "output_file_name", shortName = "outFile", doc = "Output aligned methy value files name", required = false)
    public String gchFile = null;
	
	@Output(fullName = "datapoint_file_name", shortName = "dataPoint", doc = "Output cytosine data point files name", required = false)
    public String gchDataPoint = null;
	
//	@Argument(fullName = "feature_name", shortName = "feature", doc = "Feature name provide in -B:<name>,<type> <filename> option", required = false)
 //   public String feature = null;
	
	@Argument(fullName = "search_distance_to_feature", shortName = "distance", doc = "define the distance before or after feature, default: 2000", required = false)
    public int distance = 2000;
	
	@Argument(fullName = "bin_size", shortName = "binSize", doc = "define the bin size when sliding window. default: 1", required = false)
    public int binSize = 1;
	
	@Argument(fullName = "minium_CT_reads_count", shortName = "minCTdepth", doc = "minium number of CT reads should contained to calculate methylation value, default: 1", required = false)
    public int minCTdepth = 1;
	
	//@Argument(fullName = "space_before_feature", shortName = "before", doc = "define the space before feature to detect", required = false)
  //  public int before = 2000;
	
	///@Argument(fullName = "space_after_feature", shortName = "after", doc = "define the space after feature to detect", required = false)
   // public int after = 2000;
	
	@Argument(fullName = "enable_orientation", shortName = "orientated", doc = "orientated by strand or not, default: not orientated", required = false)
    public boolean orientated = false;
	
	@Argument(fullName = "avoid_overlap_itself", shortName = "noItself", doc = "when aligned feature are overlapped with input value, avoid the alignment, default: not enabled", required = false)
    public boolean noItself = false;
	
	@Argument(fullName = "bed_format", shortName = "bedFormat", doc = "define the bed format of input methylation value: " +
			"[BED6PLUS2: TCGA bed 6+2 format, whose score column is 0-1000 to represent methylation value, this is default option;" +
			"BED3PLUS2: bed 3+2 format, whose score column is 0-100 to represent methylation value;" +
			"ENCODE: Encode bed 9 format, whose score column is 0-1000 to represent methylation value]", required = false)
    public BedFormat bedFormat = BedFormat.BED6PLUS2;
	
	@Argument(fullName = "motif_alignment_type", shortName = "alignmentType", doc = "motif aligned at FiveEnd, ThreeEnd or Center, default: align to center", required = false)
    public MotifAlignmentType alignmentType = MotifAlignmentType.Center;
	
	//output three files matrix file, contained GCH, WCG and HCG
	private bedObjectWriterImp gchWriter = null;

	
	private bedObjectWriterImp gchDataWriter = null;
	
	private boolean inFeature = false;
	private boolean writtenObject = false;
	
	private LinkedList<Double> tmpMethyValueListGch = null;
	
	private LinkedList<Integer> tmpDataPointListGch = null;
	
	private ReferenceOrderedDataSource rodIt = null;
	
	
	private String chr = null;
	private int bedStart = 0;
	private int bedEnd = 0;
	private SimpleBEDFeature bed = null;
	private Strand strand = Strand.NONE;

	
	
	
	
	/**
	 * 
	 */
	public void initialize(){
		 File fn1 = null;
		 if(gchFile != null)
				fn1 = new File(gchFile);
		 File fn4= null;
		 if(gchDataPoint != null)
			fn4 = new File(gchDataPoint);
		 if(gchFile != null)
			 gchWriter = new bedObjectWriterImp(fn1);

		 if(gchDataPoint != null)
			 gchDataWriter = new bedObjectWriterImp(fn4);
		 
		 tmpMethyValueListGch = new LinkedList<Double>();

		 tmpDataPointListGch = new LinkedList<Integer>();
		 rodIt = getToolkit().getRodDataSources().get(1);
	}

	

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.LocusWalker#map(org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker, org.broadinstitute.sting.gatk.contexts.ReferenceContext, org.broadinstitute.sting.gatk.contexts.AlignmentContext)
	 */
	@Override
	public Boolean map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
	//	if(tracker.getValues(values).isEmpty()){
	//		return null;
	//	}
	//	else{
	//		System.err.println(tracker.getValues(values).get(0).getStart() + "\t" + tracker.getValues(values).get(0).getScore() + "\t" + tracker.getValues(values).get(0).getStrand());
	//		System.err.println(tracker.getValues(values).get(0).getDescription() + "\t" + tracker.getValues(values).get(0).getLink() + "\t" + tracker.getValues(values).get(0).getName());
	//		System.err.println(tracker.getValues(values).get(0).getType() + "\t" + tracker.getValues(values).get(0).getColor() + "\t" + tracker.getValues(values).get(0).getExons());
	//	}
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
   					int featureAlignStart = (bedTmp.getStart() + bedTmp.getEnd())/2;
					 
					 if(alignmentType == MotifAlignmentType.FiveEnd){
						 featureAlignStart = bedTmp.getStart();
						 if(bedTmp.getStrand() == Strand.NEGATIVE){
							 featureAlignStart = bedTmp.getEnd();
						 }
					 }
					 else if(alignmentType == MotifAlignmentType.ThreeEnd){
						 featureAlignStart = bedTmp.getEnd();
						 if(bedTmp.getStrand() == Strand.NEGATIVE){
							 featureAlignStart = bedTmp.getStart();
						 }
					 }
					 
       				 if(loc.distance(getToolkit().getGenomeLocParser().createGenomeLoc(bedTmp.getChr(), featureAlignStart, featureAlignStart)) <= distance){
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
   				int featureAlignStart = (bed.getStart() + bed.getEnd())/2;
			 
   				if(alignmentType == MotifAlignmentType.FiveEnd){
   					featureAlignStart = bed.getStart();
   					if(bed.getStrand() == Strand.NEGATIVE){
   						featureAlignStart = bed.getEnd();
   					}
   				}
   				else if(alignmentType == MotifAlignmentType.ThreeEnd){
   					featureAlignStart = bed.getEnd();
   					if(bed.getStrand() == Strand.NEGATIVE){
   						featureAlignStart = bed.getStart();
   					}
   				}	 
			 
   					if(bed != null && loc.distance(getToolkit().getGenomeLocParser().createGenomeLoc(bed.getChr(), featureAlignStart, featureAlignStart)) > distance){

    	    		 inFeature = false;
    	    		 writtenObject = false;
    	    		// System.err.println(bed.getStart());
    	    		// break;
    				 }
    		}
   			 
           
   	 }
   	 
   	
   	 rodIt.close(locRodIt);
   	
	     if(!inFeature){
	    	 if(!writtenObject && (!tmpMethyValueListGch.isEmpty() || !tmpDataPointListGch.isEmpty())){
	    		// System.err.println(tmpMethyValueListGch.size() + "\t" + tmpMethyValueListWcg.size());
	    		 if(tmpMethyValueListGch.size() == distance * 2 + 1 || tmpDataPointListGch.size() == distance * 2 + 1){
	    			 if(gchFile != null){
	    				 bedObject bedLineGch = new bedObject(chr, bedStart, bedEnd, strand, (List)tmpMethyValueListGch); //chr, start, end, strand, aveMethyNDR, gchNumNDR, gchDepthNDR, gchCTdepthNDR, aveMethyLinker, gchNumLinker, gchDepthLinker, gchCTdepthLinker
	 		    		// System.err.println("size:" + tmpMethyValueListGch.size());
	 	    			 gchWriter.add(bedLineGch);
	    			 }
	    			
		    		 if(gchDataPoint != null){
		    			 bedObject bedLineGch2 = new bedObject(chr, bedStart, bedEnd, strand, (List)tmpDataPointListGch);
			    		 gchDataWriter.add(bedLineGch2);
		    		 }
		    		 
	    		 }
	    		 if(gchFile != null)
	    			 tmpMethyValueListGch.clear();
	    		
	    		 if(gchDataPoint != null)
	    			 tmpDataPointListGch.clear();
	    		 
	    		 writtenObject = true;
	    		 logger.info(chr + "\t" + bedStart + "\t" + bedEnd);
	    	 }
	    	
	    	 return null;
	     }
	     else{

	    	 List<BEDFeature> bedValues = tracker.getValues(values);
	    	 double methy = Double.NaN;
	    	 int cov;
	    	 if(bedValues.isEmpty()){
	    		 methy = Double.NaN;
	    		 cov = -1;
	    	 }
	    	 else{
	    		 if(bedFormat == BedFormat.BED6PLUS2){
	    			 methy = (double)bedValues.get(0).getScore()/(double)1000;
	    		 }
	    		 else if(bedFormat == BedFormat.ENCODE){
	    			 methy = (double)bedValues.get(0).getScore()/(double)1000;
	    		 }
	    		 else if(bedFormat == BedFormat.BED3PLUS2){
	    			 methy = (double)bedValues.get(0).getScore()/(double)100;
	    		 }
	    		 else{
	    			 try {
						throw new Exception("not recognized bed file format");
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
	    		 }
	    		 if(noItself && (bedValues.get(0).getChr().equalsIgnoreCase(chr) && bedValues.get(0).getStart()==bedStart)){ // bug? should be bedStart-1?
	    			 cov = -1;
	    		 }
	    		 cov = 1;
	    	 }
	    	 if(gchFile != null)
	    		 addContextToList(methy, strand, tmpMethyValueListGch);
	    	 if(gchDataPoint != null)
	    		 addCoverageToList(cov, strand, tmpDataPointListGch);
	     }
		return null;

		
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduceInit()
	 */
	@Override
	public Boolean reduceInit() {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Boolean reduce(Boolean value, Boolean sum) {
		// TODO Auto-generated method stub
		return null;
	}
	
	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.TreeReducible#treeReduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Boolean treeReduce(Boolean lhs, Boolean rhs) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public void onTraversalDone(Datum result) {
		if(gchFile != null)
			gchWriter.close();
		
   	 	if(gchDataPoint != null)
   	 		gchDataWriter.close();
   	 	
		logger.info("Finished!");
	}
	
	private void addContextToList(double methy, Strand strand, LinkedList<Double> list){
		if(orientated){
			if(strand == Strand.NEGATIVE){
				
					list.offerFirst(methy);
			}
			else{
				
					list.offerLast(methy);
			}
			
		}
		else{
				list.offerLast(methy);
		}
	}
	
	private void addCoverageToList(int coverage, Strand strand, LinkedList<Integer> list){
		if(orientated){
			if(strand == Strand.NEGATIVE){
				if(coverage == -1){
					list.offerFirst(0);
				}
				else{
					list.offerFirst(coverage);
				}
			}
			else{
				if(coverage == -1){
					list.offerLast(0);
				}
				else{
					list.offerLast(coverage);
				}
			}
			
		}
		else{
			if(coverage == -1){
				list.offerLast(0);
			}
			else{
				list.offerLast(coverage);
			}
		}
	}
	
	private enum MotifAlignmentType{
		FiveEnd,
		ThreeEnd,
		Center
	}
	
	private enum BedFormat{
		ENCODE, //score is 0-100%, but a little different..
		BED6PLUS2, // score is 0-1000%
		BED3PLUS2 // score is 0-100%
	}

}
