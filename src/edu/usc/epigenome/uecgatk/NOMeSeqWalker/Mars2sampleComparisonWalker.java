/**
 * output MARs with status changed marked:
 * 
 * NoChange: Not changed
 * NewMARs: new MARs in sample 2
 * LostMARs: MARs lost in sample 2
 * LenIncr: MARs increase in sample 2
 * LenDecr: MARs decrease in sample 2
 * BoundShift:MARs boundary shift
 * Vague: changes not clear in sample 2
 * SPR: New SPR in sample2
 * MonoNuc: New MPR in sample2
 * MPR: New MPR in sample2
 * SPR_lost: SPR lost in sample2
 * MonoNuc_lost: MonoNuc lost in sample2
 * MPR_lost:MPR lost in sample2
 * 
 */
package edu.usc.epigenome.uecgatk.NOMeSeqWalker;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.TreeSet;
import java.util.List;

import org.broad.tribble.bed.BEDFeature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.PartitionBy;
import org.broadinstitute.sting.gatk.walkers.PartitionType;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import edu.usc.epigenome.uecgatk.YapingWriter.bedObject;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObjectWriterImp;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Feb 21, 2013 12:29:28 PM
 * 
 */
@By(DataSource.REFERENCE)
@Requires({DataSource.REFERENCE, DataSource.REFERENCE_BASES})
@PartitionBy(PartitionType.LOCUS)
public class Mars2sampleComparisonWalker extends LocusWalker<Long, Long> implements TreeReducible<Long> {
	/**
	 * -L will be sample A and sample B's union segment (significant MARs' union)
	 */
	@Input(fullName="mar_bed_1", shortName = "mar1", doc="significant MAR bed files in sample 1", required=true) // this is just use MARs as an example, if use MPR, then it is the same case..
	 public RodBinding<BEDFeature> mar1;
	
//	@Input(fullName="mar_seg_1", shortName = "marSeg1", doc="All MAR bed files in sample 1", required=true) 
//	 public RodBinding<BEDFeature> marSeg1;
	
	@Input(fullName="hmm_state_1", shortName = "hmm1", doc="HMM states in each GCH in sample 1", required=false) 
	 public RodBinding<BEDFeature> hmm1;
	
	@Input(fullName="mar_bed_2", shortName = "mar2", doc="significant MAR bed files in sample 2", required=false)
	 public RodBinding<BEDFeature> mar2;
	
	@Input(fullName="mar_seg_2", shortName = "marSeg2", doc="All MAR bed files in sample 2", required=false) 
	 public RodBinding<BEDFeature> marSeg2;
	
	@Input(fullName="hmm_state_2", shortName = "hmm2", doc="HMM states in each GCH in sample 2", required=false) 
	 public RodBinding<BEDFeature> hmm2;

	@Input(fullName="mpr_bed_1", shortName = "mpr1", doc="significant MPR bed files in sample 1", required=false)
	 public RodBinding<BEDFeature> mpr1;

	@Input(fullName="mpr_bed_2", shortName = "mpr2", doc="significant MPR bed files in sample 2", required=false)
	 public RodBinding<BEDFeature> mpr2;

	
	@Input(fullName="mpr_seg_2", shortName = "mprSeg2", doc="All MPR bed files in sample 2", required=false) 
	 public RodBinding<BEDFeature> mprSeg2;

	
	@Argument(fullName="min_data_point", shortName = "minCT", doc="minimum number of CT reads in the segment to be considered as a callable region. Default: 5", required=false)
	public int minCT = 5;
	
	@Argument(fullName="min_Gch", shortName = "minGch", doc="minimum number of confident GCHs in the bed_B segment to be considered as a callable region. Default: 2", required=false)
	public int minGch = 2;
	
	@Argument(fullName="mar_size_change", shortName = "msc", doc="minimum fraction of mar's length changed are considered as MAR size change. Default: 0.5", required=false)
	public double msc = 0.5;

	@Output(doc="Output MARs change status segment bed file", required=true) 
	public String outFile;
	
	private segStatus stat_previous_point = segStatus.NotKnown;
	
	private GenomeLocParser gParser = null;
	
	private bedObjectWriterImp writer = null;

	private Datum segment = null;
	
	public void initialize() {
		//rodIt = getToolkit().getRodDataSources().get(0);
		gParser = new GenomeLocParser(getToolkit().getMasterSequenceDictionary());
		File fn1 = new File(outFile);
		writer = new bedObjectWriterImp(fn1);
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.LocusWalker#map(org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker, org.broadinstitute.sting.gatk.contexts.ReferenceContext, org.broadinstitute.sting.gatk.contexts.AlignmentContext)
	 */
	@Override
	public Long map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
		if (tracker == null)
			return null;
		List<BEDFeature> mar1Values = tracker.getValues(mar1);
		//List<BEDFeature> marSeg1Values = tracker.getValues(marSeg1);
		List<BEDFeature> mpr1Values = tracker.getValues(mpr1);
		List<BEDFeature> mar2Values = tracker.getValues(mar2);
		//List<BEDFeature> marSeg2Values = tracker.getValues(marSeg2);
		List<BEDFeature> mpr2Values = tracker.getValues(mpr2);
		//List<BEDFeature> mprSeg2Values = tracker.getValues(mprSeg2);
		
		List<BEDFeature> hmm1Values = tracker.getValues(hmm1);
		List<BEDFeature> hmm2Values = tracker.getValues(hmm2);
		/*
		if(!mar1Values.isEmpty()){
			System.err.println("mar1");
			System.err.println(mar1Values.get(0).toString());
			 System.err.println(mar1Values.get(0).getChr());
			 System.err.println(mar1Values.get(0).getStart());
			 System.err.println(mar1Values.get(0).getEnd());
		}
		if(!mar2Values.isEmpty()){
			System.err.println("mar2");
			System.err.println(mar2Values.get(0).toString());
			 System.err.println(mar2Values.get(0).getChr());
			 System.err.println(mar2Values.get(0).getStart());
			 System.err.println(mar2Values.get(0).getEnd());
		}
		if(!mpr1Values.isEmpty()){
			System.err.println("mpr1");
			System.err.println(mpr1Values.get(0).toString());
			 System.err.println(mpr1Values.get(0).getChr());
			 System.err.println(mpr1Values.get(0).getStart());
			 System.err.println(mpr1Values.get(0).getEnd());
		}
		if(!mpr2Values.isEmpty()){
			System.err.println("mpr2");
			System.err.println(mpr2Values.get(0).toString());
			 System.err.println(mpr2Values.get(0).getChr());
			 System.err.println(mpr2Values.get(0).getStart());
			 System.err.println(mpr2Values.get(0).getEnd());
		}
		*/
		
		if(!mar1Values.isEmpty() || !mar2Values.isEmpty()){
			
			if(segment == null){
				segment = new Datum(minGch, minCT);
			}
			if(!mar1Values.isEmpty())
				segment.add(new Seg(mar1Values.get(0)),segment.segInSample1);
			if(!mar2Values.isEmpty())
				segment.add(new Seg(mar2Values.get(0)),segment.segInSample2);
			if(!hmm1Values.isEmpty())
				segment.add(new HmmState(hmm1Values.get(0)),segment.stateInSample1);
			if(!hmm2Values.isEmpty())
				segment.add(new HmmState(hmm2Values.get(0)),segment.stateInSample2);
			if(!mpr1Values.isEmpty())
				segment.add(new Seg(mpr1Values.get(0)),segment.mprInSample1);
			if(!mpr2Values.isEmpty())
				segment.add(new Seg(mpr2Values.get(0)),segment.mprInSample2);

				
		}
		else{
			if(segment != null){
				segment.merge();
				segment.segChangeIdentify(writer, segment.segMergedInSample1, segment.segMergedInSample2, segment.stateInSample1, segment.stateInSample2, segment.mprInSample1, segment.mprInSample2);
				segment.clear();
				segment = null;
			}
		}
		return null;
	}
	
	 public void onTraversalDone(Long sum) {
		 
		 writer.close();
		 logger.info("Finished!");
		 
		 
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

	 private enum segStatus{
		 MARs,
		 MARs_seg,
		 MPRs,
		 MPRs_seg,
		 No_overlap,
		 NotKnown
	 }

	 private class Datum{
		 TreeSet<Seg> segInSample1;
		 TreeSet<Seg> segInSample2;
		 TreeSet<HmmState> stateInSample1;
		 TreeSet<HmmState> stateInSample2;
		 TreeSet<Seg> segMergedInSample1;
		 TreeSet<Seg> segMergedInSample2;
		 
		 TreeSet<Seg> mprInSample1;
		 TreeSet<Seg> mprInSample2;

		 
		 int minCT;
		 int minGch;
		 
		 public Datum(int minCT, int minGch){ // this is the original MARs segment in sample 1 for comparison
			 segInSample1 = new TreeSet<Seg>();
			 segInSample2 = new TreeSet<Seg>();
			 stateInSample1 = new TreeSet<HmmState>();
			 stateInSample2 = new TreeSet<HmmState>();
			 
			 segMergedInSample1 = new TreeSet<Seg>();
			 segMergedInSample2 = new TreeSet<Seg>();
			 
			 mprInSample1 = new TreeSet<Seg>();
			 mprInSample2 = new TreeSet<Seg>();

			 this.minCT = minCT;
			 this.minGch = minGch;
		 }
		 
		 public void add(Seg data, TreeSet<Seg> segInSample){
			 segInSample.add(data); 
		 }
		 
		 public void add(HmmState data, TreeSet<HmmState> stateInSample){
			 stateInSample.add(data); 

		 }
		 
		 public int numGchInInterval(GenomeLoc loc, TreeSet<HmmState> stateInSample){
			int numGch = 0;
			for(HmmState state: stateInSample){
				if(loc.containsP(state.position)){
					numGch++;
				}
			}
			 return numGch;

		 }
		 
		 public int numCTInInterval(GenomeLoc loc, TreeSet<HmmState> stateInSample){
			 int numCT = 0;
				for(HmmState state: stateInSample){
					if(loc.containsP(state.position)){
						numCT = numCT + state.numC + state.numT;
					}
				}
				 return numCT;
		 }
		 
		 public boolean enoughDataInInterval(GenomeLoc loc, TreeSet<HmmState> stateInSample){
			 int numCT = numCTInInterval(loc, stateInSample);
			int numGch = numGchInInterval(loc, stateInSample);
			return 	(numCT >= minCT && numGch >= minGch);
			
		 }
		 
		 public boolean enoughDataInIntervalList(TreeSet<GenomeLoc> locs, TreeSet<HmmState> stateInSample){
			 int numCT = 0;
			int numGch = 0;
			for(GenomeLoc loc : locs){
				numCT += numCTInInterval(loc, stateInSample);
				numGch += numGchInInterval(loc, stateInSample);
			}
			return 	(numCT >= minCT && numGch >= minGch);
			
		 }
		 
		 public boolean overlap(GenomeLoc loc, TreeSet<Seg> segInSample){
			 for(Seg seg : segInSample){
				 if(loc.overlapsP(seg.position))
					 return true;
			 }
			 return false;

			
		 }
		 
		 public int sizeOfNotOverlap(GenomeLoc loc1, GenomeLoc loc2){
			 int size = 0;
			 if(loc1.overlapsP(loc2)){
				 if(loc1.containsP(loc2)){
					 for(GenomeLoc loc : loc1.subtract(loc2)){
						 size += loc.size();
					 }
				 }
				 else if(loc2.containsP(loc1)){
					 for(GenomeLoc loc : loc2.subtract(loc1)){
						 size += loc.size();
					 }
				 }
				 else{
					 for(GenomeLoc loc : loc1.subtract(loc2)){
						 size += loc.size();
					 }
					 for(GenomeLoc loc : loc2.subtract(loc1)){
						 size += loc.size();
					 }
				 }
			 }
			 else{
				 size = loc1.size() + loc2.size();
			 }

			 return size;
		 }
		
		 
		 public void clear(){
			 segInSample1.clear();
			 segInSample2.clear();
			 stateInSample1.clear();
			 stateInSample2.clear();
			 segMergedInSample1.clear();
			 segMergedInSample2.clear();
		 }
		 
		 public void merge(){
			 if(segInSample1.size() > 1){
				 segMergedInSample1 = merge(segInSample1, stateInSample1);
			 }
			 else{
				 segMergedInSample1.addAll(segInSample1);
			 }
				 
			 if(segInSample2.size() > 1){
				 segMergedInSample2 = merge(segInSample2, stateInSample2);
			 }
			 else{
				 segMergedInSample2.addAll(segInSample2); 
			 }
				 
		 }
		 
		 public TreeSet<Seg> merge(TreeSet<Seg> segInSample, TreeSet<HmmState> stateInSample){
			
			 
			 Iterator<Seg> it = segInSample.iterator();
			Seg preSeg = null;
			GenomeLoc preLoc = null;
			TreeSet<Seg> segMerged = new TreeSet<Seg>();
			//if(segInSample.size() <=1){
			//	segMerged.addAll(segInSample);
			//	return segMerged;
			//}
				
			while(it.hasNext()){
				if(preLoc == null){
					preSeg = it.next();
					preLoc = preSeg.position;
				}
				else{
					Seg currentSeg = it.next();
					GenomeLoc interLoc = gParser.createGenomeLoc(preLoc.getContig(),preLoc.getStop(),currentSeg.position.getStart());

					if(!enoughDataInInterval(interLoc, stateInSample)){
						Seg tmpSeg = merge(preSeg, currentSeg);
						preLoc = tmpSeg.position;
						preSeg = tmpSeg;
					}
					else{
						segMerged.add(preSeg);
						preLoc = currentSeg.position;
						preSeg = currentSeg;
					}
				}
			}
			
			 return segMerged;
		 }
		 
		 public Seg merge(Seg seg1, Seg seg2){//here, seg2 should be after seg1 in the genome location
			 String chr = seg1.position.getContig();
			 int start = Math.min(seg1.position.getStart(),seg2.position.getStart());
			 int end = Math.min(seg1.position.getStop(),seg2.position.getStop());
			 return new Seg(gParser.createGenomeLoc(chr,start,end));
			 //return new Seg(gParser.createGenomeLoc(chr,start,end),(seg1.numGch+seg2.numGch), (seg1.numC+seg2.numC), (seg1.numT+seg2.numT), Math.min(seg1.fdr, seg2.fdr), seg1.alpha1, seg2.alpha2);
		 }
/*
 * NoChange: Not changed
 * NewMARs: new MARs in sample 2
 * LostMARs: MARs lost in sample 2
 * LenIncr: MARs increase in sample 2
 * LenDecr: MARs decrease in sample 2
 * BoundShift:MARs boundary shift
 * Vague: changes not clear in sample 2
 * SPR: New SPR in sample2
 * MonoNuc: New MPR in sample2
 * MPR: New MPR in sample2
 * SPR_lost: SPR lost in sample2
 * MonoNuc_lost: MonoNuc lost in sample2
 * MPR_lost:MPR lost in sample2
 */

		 /*5th: overlap or not
		 * 6th: how changes detail
		 * 7th: overlapped MARs location
		 * 8th: if there are MPR in another sample, how is the location
		 * 9th: how many GCH in  sample1's not overlapped region
		 * 10th: how many numCT in  sample1's not overlapped region
		 * 11th: how many GCH in  sample2's not overlapped region
		 * 12th: how many numCT in  sample2's not overlapped region
		 */
		 public void segChangeIdentify(bedObjectWriterImp writer, TreeSet<Seg> segMergedInSample1, TreeSet<Seg> segMergedInSample2, TreeSet<HmmState> stateInSample1, TreeSet<HmmState> stateInSample2, TreeSet<Seg> mprInSample1, TreeSet<Seg> mprInSample2){
			//System.err.println(segMergedInSample1.size() + "\t" + segMergedInSample2.size());
			 if(segMergedInSample1.size() == 1 && segMergedInSample2.size() == 0){ //LostMARs, MAR in sample 1 without enough data points will not be in output
				//if(enoughDataInInterval(segMergedInSample1.first().position, stateInSample2)){ //if incorporate BetaBB
				if(enoughDataInInterval(segMergedInSample1.first().position, stateInSample2) && enoughDataInInterval(segMergedInSample1.first().position, stateInSample1)){ //if incorporate BetaBB
					List<Object> content = new ArrayList<Object>();
					content.add("NoOverlap");
					if(mprInSample2.size()!=0){
						
						content.add("LostMAR;NewMPR");
						String mprLocs = "";
						for(Seg mpr : mprInSample2){
							mprLocs += (mpr.position.toString() + ";");
						}	
						content.add(mprLocs);
					}
					else{
						content.add("LostMAR");
						content.add("NA");
					}
					content.add(numGchInInterval(segMergedInSample1.first().position, stateInSample1));
					content.add(numCTInInterval(segMergedInSample1.first().position, stateInSample1));
					content.add(numGchInInterval(segMergedInSample1.first().position, stateInSample2));
					content.add(numCTInInterval(segMergedInSample1.first().position, stateInSample2));
					bedObject tmpLine = new bedObject(segMergedInSample1.first().position.getContig(),segMergedInSample1.first().position.getStart(),segMergedInSample1.first().position.getStop(),content);
					writer.add(tmpLine);
				}
			}
			else if(segMergedInSample2.size() == 1 && segMergedInSample1.size() == 0){ //NewMARs, MAR in sample 2 without enough data points will not be in output
				//if(enoughDataInInterval(segMergedInSample2.first().position, stateInSample1)){
				if(enoughDataInInterval(segMergedInSample2.first().position, stateInSample1) && enoughDataInInterval(segMergedInSample2.first().position, stateInSample2)){
					List<Object> content = new ArrayList<Object>();
					content.add("NoOverlap");
					if(mprInSample1.size()!=0){
						
						content.add("LostMPR;NewMAR");
						String mprLocs = "";
						for(Seg mpr : mprInSample1){
							mprLocs += (mpr.position.toString() + ";");
						}
						content.add(mprLocs);
					}
					else{
						content.add("NewMAR");
						content.add("NA");
					}
					content.add(numGchInInterval(segMergedInSample2.first().position, stateInSample1));
					content.add(numCTInInterval(segMergedInSample2.first().position, stateInSample1));
					content.add(numGchInInterval(segMergedInSample2.first().position, stateInSample2));
					content.add(numCTInInterval(segMergedInSample2.first().position, stateInSample2));
					bedObject tmpLine = new bedObject(segMergedInSample2.first().position.getContig(),segMergedInSample2.first().position.getStart(),segMergedInSample2.first().position.getStop(),content);
					writer.add(tmpLine);
				}
			}
			else if(segMergedInSample2.size() != 0 && segMergedInSample1.size() != 0){ //MARsOverlap
				if(segMergedInSample1.size() == 1 && segMergedInSample2.size() == 1){ // 1 to 1
					GenomeLoc locInSample1 = segMergedInSample1.first().position;
					GenomeLoc locInSample2 = segMergedInSample2.first().position;
					List<Object> content = new ArrayList<Object>();
					content.add("Overlap");
					if(locInSample1.equals(locInSample2)){ //no change at all

						content.add("MARnoChange");
						content.add("NA");
						content.add(0);
						content.add(0);
						content.add(0);
						content.add(0);
						bedObject tmpLine = new bedObject(segMergedInSample1.first().position.getContig(),segMergedInSample1.first().position.getStart(),segMergedInSample1.first().position.getStop(),content);
						writer.add(tmpLine);
					}
					else if(locInSample1.containsP(locInSample2)){//Length decrease
						
						if((double)sizeOfNotOverlap(locInSample1, locInSample2)/(double)locInSample1.size() > msc){ //Length decrease dramatically
							int trueOverlapSize = 0;
							boolean mprFlag = false;
							int numGchInSample1 = 0;
							int numCTInSample1 = 0;
							int numGchInSample2 = 0;
							int numCTInSample2 = 0;
							for(GenomeLoc loc : locInSample1.subtract(locInSample2)){
								numGchInSample1 += numGchInInterval(loc, stateInSample1);
								numGchInSample1 += numCTInInterval(loc, stateInSample1);
								numGchInSample2 += numGchInInterval(loc, stateInSample2);
								numGchInSample2 += numCTInInterval(loc, stateInSample2);
								if(enoughDataInInterval(loc, stateInSample2)){ //Length decrease due to lack of data point in sample 2?
									trueOverlapSize += loc.sizeOfOverlap(locInSample1);
									
										if(overlap(loc, mprInSample2)){//Length decrease due to new MPR in sample 2?
											mprFlag = true;
										}
								}
							}
							if((double)trueOverlapSize/(double)locInSample1.size() > msc){
								if(mprFlag){
									content.add("MARlenDecrease;NewMPR");
									String mprLocs = "";
									for(Seg mpr : mprInSample2){
										mprLocs += (mpr.position.toString() + ";");
									}
									content.add(mprLocs);
									
								}
								else{
									content.add("MARlenDecrease");
									content.add("NA");
								}
								content.add(numGchInSample1);
								content.add(numCTInSample1);
								content.add(numGchInSample2);
								content.add(numCTInSample2);
								bedObject tmpLine = new bedObject(segMergedInSample1.first().position.getContig(),segMergedInSample1.first().position.getStart(),segMergedInSample1.first().position.getStop(),content);
								writer.add(tmpLine);
							}
							else{// length decrease because of lack of data point in the boundary
								content.add("LackData");
								content.add("NA");
								content.add(numGchInSample1);
								content.add(numCTInSample1);
								content.add(numGchInSample2);
								content.add(numCTInSample2);
								bedObject tmpLine = new bedObject(segMergedInSample1.first().position.getContig(),segMergedInSample1.first().position.getStart(),segMergedInSample1.first().position.getStop(),content);
								writer.add(tmpLine);
							}
						}
						else{//Length decrease NOT dramatically
							content.add("MARnoChange");
							content.add("NA");
							content.add(0);
							content.add(0);
							content.add(0);
							content.add(0);
							bedObject tmpLine = new bedObject(segMergedInSample1.first().position.getContig(),segMergedInSample1.first().position.getStart(),segMergedInSample1.first().position.getStop(),content);
							writer.add(tmpLine);
						}
						

					}
					else if(locInSample2.containsP(locInSample1)){//Length increase
						if((double)sizeOfNotOverlap(locInSample1, locInSample2)/(double)locInSample1.size() > msc){ //Length increase dramatically
							int trueOverlapSize = 0;
							boolean mprFlag = false;
							int numGchInSample1 = 0;
							int numCTInSample1 = 0;
							int numGchInSample2 = 0;
							int numCTInSample2 = 0;
							for(GenomeLoc loc : locInSample2.subtract(locInSample1)){
								numGchInSample1 += numGchInInterval(loc, stateInSample1);
								numGchInSample1 += numCTInInterval(loc, stateInSample1);
								numGchInSample2 += numGchInInterval(loc, stateInSample2);
								numGchInSample2 += numCTInInterval(loc, stateInSample2);
								if(enoughDataInInterval(loc, stateInSample1)){ //Length increase due to lack of data point in sample 1?
									trueOverlapSize += loc.sizeOfOverlap(locInSample2);
									
									
										if(overlap(loc, mprInSample1)){//Length increase due to MPR lost in sample 1?
											mprFlag = true;
										}
								}
							}
							if((double)trueOverlapSize/(double)locInSample1.size() > msc){
								if(mprFlag){
									content.add("MARlenIncrease;LostMPR");
									String mprLocs = "";
									for(Seg mpr : mprInSample1){
										mprLocs += (mpr.position.toString() + ";");
									}
									content.add(mprLocs);
								}
								else{
									content.add("MARlenIncrease");
									content.add("NA");
								}
								content.add(numGchInSample1);
								content.add(numCTInSample1);
								content.add(numGchInSample2);
								content.add(numCTInSample2);
								bedObject tmpLine = new bedObject(segMergedInSample1.first().position.getContig(),segMergedInSample1.first().position.getStart(),segMergedInSample1.first().position.getStop(),content);
								writer.add(tmpLine);
							}
							else{// length increase because of lack of data point in the boundary
								content.add("LackData");
								content.add("NA");
								content.add(numGchInSample1);
								content.add(numCTInSample1);
								content.add(numGchInSample2);
								content.add(numCTInSample2);
								bedObject tmpLine = new bedObject(segMergedInSample1.first().position.getContig(),segMergedInSample1.first().position.getStart(),segMergedInSample1.first().position.getStop(),content);
								writer.add(tmpLine);
							}
						}
						else{//Length increase NOT dramatically
							content.add("MARnoChange");
							content.add("NA");
							content.add(0);
							content.add(0);
							content.add(0);
							content.add(0);
							bedObject tmpLine = new bedObject(segMergedInSample1.first().position.getContig(),segMergedInSample1.first().position.getStart(),segMergedInSample1.first().position.getStop(),content);
							writer.add(tmpLine);
						}
					}
					else if(locInSample2.overlapsP(locInSample1)){//partially overlapped, it is defined as boundary shift, when the center of sample 2 is out of sample 1' region, vice versa.
						GenomeLoc locUniqInSample1 = locInSample1.subtract(locInSample2).get(0);
						GenomeLoc locUniqInSample2 = locInSample2.subtract(locInSample1).get(0);
						int numGchInSample1 = numGchInInterval(locUniqInSample1,stateInSample1) + numGchInInterval(locUniqInSample2,stateInSample1);
						int numCTInSample1 = numCTInInterval(locUniqInSample1, stateInSample1) + numCTInInterval(locUniqInSample2, stateInSample1);
						int numGchInSample2 = numGchInInterval(locUniqInSample1,stateInSample2) + numGchInInterval(locUniqInSample2,stateInSample2);
						int numCTInSample2 = numCTInInterval(locUniqInSample1, stateInSample2) + numCTInInterval(locUniqInSample2, stateInSample2);
						
						if((double)locUniqInSample1.size()/(double)locInSample1.size() > msc && (double)locUniqInSample2.size()/(double)locInSample2.size() > msc){ //Boundary shift dramatically
							if(enoughDataInInterval(locUniqInSample1, stateInSample2) && enoughDataInInterval(locUniqInSample2, stateInSample1)){ //boundary shift is real
								if(overlap(locUniqInSample1,mprInSample2)){
									content.add("MARboundShift;NewMPR");
									String mprLocs = "";
									for(Seg mpr : mprInSample2){
										mprLocs += (mpr.position.toString() + ";");
									}
									content.add(mprLocs);
								}
								else if(overlap(locUniqInSample2,mprInSample1)){
									content.add("MARboundShift;LostMPR");
									String mprLocs = "";
									for(Seg mpr : mprInSample1){
										mprLocs += (mpr.position.toString() + ";");
									}
									content.add(mprLocs);
								}
								else{
									content.add("MARboundShift");
									content.add("NA");
								}
								
								
							}
							else{
								content.add("MARnoChange");
								content.add("NA");
							}
							content.add(numGchInSample1);
							content.add(numCTInSample1);
							content.add(numGchInSample2);
							content.add(numCTInSample2);
							bedObject tmpLine = new bedObject(segMergedInSample1.first().position.getContig(),segMergedInSample1.first().position.getStart(),segMergedInSample1.first().position.getStop(),content);
							writer.add(tmpLine);
							
						}
						else{
							content.add("MARnoChange");
							content.add("NA");
							content.add(numGchInSample1);
							content.add(numCTInSample1);
							content.add(numGchInSample2);
							content.add(numCTInSample2);
							bedObject tmpLine = new bedObject(segMergedInSample1.first().position.getContig(),segMergedInSample1.first().position.getStart(),segMergedInSample1.first().position.getStop(),content);
							writer.add(tmpLine);
						}
					}
				}
				else if(segMergedInSample1.size() == 1 && segMergedInSample2.size() >= 1){ // 1 to multiple
					GenomeLoc locInSample1 = segMergedInSample1.first().position;
					List<Object> content = new ArrayList<Object>();
					content.add("Overlap");
					if(overlap(locInSample1, mprInSample2)){
						content.add("MARsplit;NewMPR");
						String mprLocs = "";
						for(Seg mpr : mprInSample2){
							mprLocs += (mpr.position.toString() + ";");
						}
						content.add(mprLocs);
					}
					else{
						content.add("MARsplit");
						content.add("NA");
					}
					int numGchInSample1 = 0;
					int numCTInSample1 = 0;
					int numGchInSample2 = 0;
					int numCTInSample2 = 0;
					for(Seg seg : segMergedInSample2){
						if(locInSample1.overlapsP(seg.position)){
							for(GenomeLoc loc : locInSample1.subtract(seg.position)){
								numGchInSample1 += numGchInInterval(loc, stateInSample1);
								numCTInSample1 += numCTInInterval(loc, stateInSample1);
								numGchInSample2 += numGchInInterval(loc, stateInSample2);
								numCTInSample2 += numCTInInterval(loc, stateInSample2);
							}
						}
						
					}
					content.add(numGchInSample1);
					content.add(numCTInSample1);
					content.add(numGchInSample2);
					content.add(numCTInSample2);
					bedObject tmpLine = new bedObject(segMergedInSample1.first().position.getContig(),segMergedInSample1.first().position.getStart(),segMergedInSample1.first().position.getStop(),content);
					writer.add(tmpLine);
				}
				else if(segMergedInSample1.size() >= 1 && segMergedInSample2.size() == 1){ // multiple to 1
					for(Seg seg1 : segMergedInSample1){ //output of each MARs in sample 1
						GenomeLoc locInSample2 = segMergedInSample2.first().position;
						List<Object> content = new ArrayList<Object>();
						content.add("Overlap");
						
						if(overlap(locInSample2, mprInSample1)){
							content.add("MARmerge;LostMPR");
							String mprLocs = "";
							for(Seg mpr : mprInSample1){
								mprLocs += (mpr.position.toString() + ";");
							}
							content.add(mprLocs);
						}
						else{
							content.add("MARmerge");
							content.add("NA");
						}
						int numGchInSample1 = 0;
						int numCTInSample1 = 0;
						int numGchInSample2 = 0;
						int numCTInSample2 = 0;
						for(Seg seg : segMergedInSample1){
							if(segMergedInSample2.first().position.overlapsP(seg.position)){
								for(GenomeLoc loc : segMergedInSample2.first().position.subtract(seg.position)){
									numGchInSample1 += numGchInInterval(loc, stateInSample1);
									numCTInSample1 += numCTInInterval(loc, stateInSample1);
									numGchInSample2 += numGchInInterval(loc, stateInSample2);
									numCTInSample2 += numCTInInterval(loc, stateInSample2);
								}
							}
							
						}
						content.add(numGchInSample1);
						content.add(numCTInSample1);
						content.add(numGchInSample2);
						content.add(numCTInSample2);
						bedObject tmpLine = new bedObject(seg1.position.getContig(),seg1.position.getStart(),seg1.position.getStop(),content);
						writer.add(tmpLine);
					}
					
				}
				else{ // multiple to multiple
					for(Seg seg1 : segMergedInSample1){
						TreeSet<Seg> newSeg1List = new TreeSet<Seg>();
						newSeg1List.add(seg1);
						segChangeIdentify(writer,newSeg1List, segMergedInSample2, stateInSample1, stateInSample2, mprInSample1, mprInSample2);
					}
				}
			}
		 }
		 
	 }
	 
	 private class Seg 
	 implements Comparable<Seg>, Comparator<Seg>{
		 public GenomeLoc position;
	//	 public int numGch;
	//	 public int numC;
	//	 public int numT;
	//	 public double fdr;
	//	 public double alpha1;
	//	 public double alpha2;
		 
		 public Seg(BEDFeature value){
			 
			 position = gParser.createGenomeLoc(value.getChr(), value.getStart(), value.getEnd()); //take care of 0-based, 1-based here..
		//	 String data = value.getDescription();
		//	 String[] datas = data.split(";");
		//	 numGch = Integer.parseInt(datas[0]);
		//	 numC = Integer.parseInt(datas[1]);
		//	 numT = Integer.parseInt(datas[2]);
		//	 fdr = Double.parseDouble(datas[3]);
		//	 alpha1 = Double.parseDouble(datas[4]);
		//	 alpha2 = Double.parseDouble(datas[5]);
		 }
		 
		 public Seg(GenomeLoc loc){
			 position = loc;
		 }
		 
		 public Seg(GenomeLoc loc, int numGch, int numC, int numT, double fdr, double alpha1, double alpha2){
			 position = loc;
		//	 this.numGch = numGch;
		//	 this.numC = numC;
		//	 this.numT = numT;
		//	 this.fdr = fdr;
		//	 this.alpha1 = alpha1;
		//	 this.alpha2 = alpha2;
		 }

		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(Seg that) {
			if(this.position.equals(that.position)){
				return 0;
			}
			else if(this.position.isBefore(that.position)){
				return -1;
			}
			else if(this.position.isPast(that.position)){
				return 1;
			}
			return 0;
		}

		/* (non-Javadoc)
		 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
		 */
		@Override
		public int compare(Seg arg0, Seg arg1) {
			if(arg0.position.equals(arg1.position)){
				return 0;
			}
			else if(arg0.position.isBefore(arg1.position)){
				return -1;
			}
			else if(arg0.position.isPast(arg1.position)){
				return 1;
			}
			return 0;
		}
		
		public boolean equals(Seg obj){
			return obj.position.equals(this.position);
			
		}
		 
	 }
	 
	 private class HmmState
	 implements Comparable<HmmState>, Comparator<HmmState>{
		 public GenomeLoc position;
		 public int state;
		 public int numC;
		 public int numT;

		 
		 public HmmState(BEDFeature value){
			 position = gParser.createGenomeLoc(value.getChr(), value.getStart(), value.getEnd()); //take care of 0-based, 1-based here..
			 String data = value.getName();
			 String[] datas = data.split(";");
			 state = Integer.parseInt(datas[0]);
			 numC = Integer.parseInt(datas[1]);
			 numT = Integer.parseInt(datas[2]) - numC;

		 }


		/* (non-Javadoc)
		 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
		 */
		@Override
		public int compare(HmmState arg0, HmmState arg1) {
			if(arg0.position.equals(arg1.position)){
				return 0;
			}
			else if(arg0.position.isBefore(arg1.position)){
				return -1;
			}
			else if(arg0.position.isPast(arg1.position)){
				return 1;
			}
			return 0;
		}


		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(HmmState that) {
			if(this.position.equals(that.position)){
				return 0;
			}
			else if(this.position.isBefore(that.position)){
				return -1;
			}
			else if(this.position.isPast(that.position)){
				return 1;
			}
			return 0;
		}
		 
	 }

}
