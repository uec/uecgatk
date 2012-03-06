/*
 * program measure NDR from the lef first position with sig value to the right first position with sig value. so NDR position would have a little shift (may be 15bp) to the left position. 
 */

package edu.usc.epigenome.uecgatk.YapingWalker;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Feb 23, 2012 10:48:18 AM
 * 
 */


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;


import jsc.independentsamples.MannWhitneyTest;
import jsc.independentsamples.SmirnovTest;
import jsc.independentsamples.TwoSampleTtest;
import jsc.tests.H1;

import net.sf.samtools.SAMRecord;

import org.apache.commons.math.distribution.BinomialDistributionImpl;
import org.apache.commons.math.util.MathUtils;
import org.apache.log4j.Logger;
import org.broad.tribble.annotation.Strand;
import org.broad.tribble.bed.SimpleBEDFeature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentReadFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;

import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisSNPUtils;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisulfiteDiploidSNPGenotypePriors;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisulfiteVariantCallContext;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObject;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObjectWriterImp;
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
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.pileup.AbstractReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import cern.jet.stat.Probability;
import edu.usc.epigenome.uecgatk.YapingWalker.NDRdetectNearLoci.windowsObject;
import edu.usc.epigenome.uecgatk.YapingWalker.NDRargumentCollection;

@BAQMode(QualityMode = BAQ.QualityMode.OVERWRITE_QUALS, ApplicationTime = BAQ.ApplicationTime.ON_INPUT)
@ReadFilters( {UnmappedReadFilter.class, BadMateFilter.class, NotPrimaryAlignmentReadFilter.class, DuplicateReadFilter.class} )
@Reference(window=@Window(start=-300,stop=300))
@By(DataSource.REFERENCE)
@Downsample(by=DownsampleType.NONE)
public class NDRdetectNearLoci extends LocusWalker<NDRCallContext,windowsObject> implements
		TreeReducible<windowsObject> {

	@ArgumentCollection private static NDRargumentCollection NAC = new NDRargumentCollection();
	
	@Argument(fullName = "feature_name", shortName = "feature", doc = "Feature name provide in -B:<name>,<type> <filename> option", required = false)
    public String feature = "tss";
	
	@Argument(fullName = "downstream_of_rod", shortName = "down", doc = "length to check in the downstream of rod, from 0 to +X bp", required = false)
    public int down = 0;
	
	@Argument(fullName = "upstream_of_rod", shortName = "up", doc = "length to check in the upstream of rod, from -X bp to 0 bp", required = false)
    public int up = 100;
    
    protected bedObjectWriterImp callableWindWriter = null;
	
	//private NDRdetectionEngine NDRD_engine = null;
	
	private windowsObject windows = null;
	
	private GenotypePriors genotypePriors;
	
	private ContextCondition summary = null;
	
	
	
	private boolean winPreLinkerflag = false;
	
	private boolean winPostLinkerflag = false;
	
	private boolean winStartflag = false;
	
	private boolean winEndflag = false;
	
	private boolean winInflag = false;
	
	private boolean winForceEndflag = false;
	
//	private boolean winPassFilterflag = true;
	
	private String winChr = null;
	
	private int winStartCor = -1;
			
	private int winEndCor = -1;
	
	private int tmpGchCTDepthWind = 0;
	
	private int tmpGchNumWind = 0;
	
	private int tmpGchDepthWind = 0;
	
	private int tmpGchDotWind = 0;
	
	private int tmpGchNumCDotWind = 0;
	
	private double tmpGchMethyWind = 0;

	private double tmpGchMethyWindLinker = 0;
	
	private double tmpWcgMethyWind = 0;
	
	private int tmpWcgNumWind = 0;
	
	//private double tmpHcgMethyWind = 0;
	
	private double sigValueMem = -1;
	
	
	private boolean outputFlag = false;
	
	private ArrayList<Double[]> gchListForVerboseMode = null;
	
	private ReferenceOrderedDataSource rodIt = null;
	
	public NDRdetectNearLoci() {
		// TODO Auto-generated constructor stub
		
	}

	 public static class ContextCondition {
		 	long nWindowsVisited = 0;

	       
		 	/** The number of windows called confidently (according to user threshold), contained enough confidantly GCH and enough seq depth and have at least one adjacent window owns enough confidantly GCH and enough seq depth*/
	        long nWindowsCallable = 0;
	        
	        /** The average sequence depth inside window */
	        double sumGchCTDepthWind = 0;
	        
	        /** The average sequence depth inside window */
	        double sumGchDepthWind = 0;
	     
	        /** The average sequence depth inside window */
	        double sumGchDotWind = 0;
	        
	        /** The number of Gch bases called confidently (according to user threshold), either ref or other */
	        long sumGchNumInWindCalledConfidently = 0;
	        
	        /** The sum of GCH methylation value of windows called confidently (according to user threshold), contained enough confidantly GCH and enough seq depth */
	        double sumGchMethyWindowsCalledConfidently = 0;
	        
	        /** The number of Gch bases called confidently (according to user threshold), either ref or other */
	        long sumHcgNumInWindCalledConfidently = 0;
	        
	        /** The sum of GCH methylation value of windows called confidently (according to user threshold), contained enough confidantly GCH and enough seq depth */
	        double sumHcgMethyWindowsCalledConfidently = 0;
	        
	        /** The number of Gch bases called confidently (according to user threshold), either ref or other */
	        long sumWcgNumInWindCalledConfidently = 0;
	        
	        /** The sum of GCH methylation value of windows called confidently (according to user threshold), contained enough confidantly GCH and enough seq depth */
	        double sumWcgMethyWindowsCalledConfidently = 0;
	        
	        /** The average sequence depth inside NDR window */
	        double sumGchCTDepthInNDRWind = 0;
	        
	        /** The average sequence depth inside NDR window */
	        double sumGchDepthInNDRWind = 0;
	        
	        /** The average sequence depth inside NDR window */
	        double sumGchDotInNDRWind = 0;
	     
	        /** The number of Gch bases called confidently (according to user threshold), either ref or other in NDR windows*/
	        long sumGchNumInNDRWind = 0;
	        
	        /** The number of NDR windows called confidently (according to user threshold), contained enough confidantly GCH and enough seq depth */
	        long nNDRWindowsCalledConfidently = 0;
	        
	        /** The sum of GCH methylation value of windows called confidently (according to user threshold), contained enough confidantly GCH and enough seq depth */
	        double sumGchMethyNDRWindowsCalledConfidently = 0;
	        
	        /** The average sequence depth inside NDR window */
	        double sumGchCTDepthInNDRWindLinker = 0;
	     
	        /** The average sequence depth inside NDR window */
	        double sumGchDepthInNDRWindLinker = 0;
	        
	        /** The average sequence depth inside NDR window */
	        double sumGchDotInNDRWindLinker = 0;
	        
	        /** The number of Gch bases called confidently (according to user threshold), either ref or other in NDR windows*/
	        long sumGchNumInNDRWindLinker = 0;
	        
	        /** The number of NDR windows called confidently (according to user threshold), contained enough confidantly GCH and enough seq depth */
	        long nNDRWindowsCalledConfidentlyLinker = 0;
	        
	        /** The sum of GCH methylation value of windows called confidently (according to user threshold), contained enough confidantly GCH and enough seq depth */
	        double sumGchMethyNDRWindowsCalledConfidentlyLinker = 0;
	        
	        
	        
	        double percentCallableWindowOfAll() { return (double)nWindowsCallable/nWindowsVisited;}
	        double percentGchMethyOfCallableWindows() { return (double)sumGchMethyWindowsCalledConfidently/(nWindowsCallable * percentGchNumOfCallableWindows());}
	        double percentHcgMethyOfCallableWindows() { return (double)sumHcgMethyWindowsCalledConfidently/(nWindowsCallable * percentHcgNumOfCallableWindows());}
	        double percentWcgMethyOfCallableWindows() { return (double)sumWcgMethyWindowsCalledConfidently/(nWindowsCallable * percentWcgNumOfCallableWindows());}
	        double percentGchCTDepthOfCallableWindows() { return (double)sumGchCTDepthWind/(nWindowsCallable * percentGchNumOfCallableWindows());}
	        double percentGchDepthOfCallableWindows() { return (double)sumGchDepthWind/(nWindowsCallable * percentGchNumOfCallableWindows());}
	        double percentGchDotOfCallableWindows() { return (double)sumGchDotWind/(nWindowsCallable);}
	        double percentGchNumOfCallableWindows() { return (double)sumGchNumInWindCalledConfidently/(nWindowsCallable);}
	        double percentHcgNumOfCallableWindows() { return (double)sumHcgNumInWindCalledConfidently/(nWindowsCallable);}
	        double percentWcgNumOfCallableWindows() { return (double)sumWcgNumInWindCalledConfidently/(nWindowsCallable);}
	        double percentNDRWindowOfCallableWindows() { return (double)nNDRWindowsCalledConfidently/nWindowsCallable;}
	        double percentGchMethyOfNDRWindowsCalledConfidently() { return (double)sumGchMethyNDRWindowsCalledConfidently/(nNDRWindowsCalledConfidently * percentGchNumOfNDRWindows());}
	        double percentGchCTDepthOfNDRWindows() { return (double)sumGchCTDepthInNDRWind/(nNDRWindowsCalledConfidently * percentGchNumOfNDRWindows());}
	        double percentGchDepthOfNDRWindows() { return (double)sumGchDepthInNDRWind/(nNDRWindowsCalledConfidently * percentGchNumOfNDRWindows());}
	        double percentGchDotOfNDRWindows() { return (double)sumGchDotInNDRWind/(nNDRWindowsCalledConfidently);}
	        double percentGchNumOfNDRWindows() { return (double)sumGchNumInNDRWind/(nNDRWindowsCalledConfidently);}
	        double percentGchMethyOfCallableWindowsLinker() { return (double)sumGchMethyNDRWindowsCalledConfidentlyLinker/(nNDRWindowsCalledConfidentlyLinker * percentGchNumOfCallableWindowsLinker());}
	        double percentGchCTDepthOfCallableWindowsLinker() { return (double)sumGchCTDepthInNDRWindLinker/(nNDRWindowsCalledConfidentlyLinker * percentGchNumOfCallableWindowsLinker());}
	        double percentGchDepthOfCallableWindowsLinker() { return (double)sumGchDepthInNDRWindLinker/(nNDRWindowsCalledConfidentlyLinker * percentGchNumOfCallableWindowsLinker());}
	        double percentGchDotOfCallableWindowsLinker() { return (double)sumGchDotInNDRWindLinker/(nNDRWindowsCalledConfidentlyLinker);}
	        double percentGchNumOfCallableWindowsLinker() { return (double)sumGchNumInNDRWindLinker/(nNDRWindowsCalledConfidentlyLinker);}
	 }
	 
	 public void initialize(){
	//	 NDRD_engine = new NDRdetectionEngine(getToolkit(), NAC, logger);
		// if(NAC.largeLinker){
		//	 modifyIntervals();
		// }
		 
		 genotypePriors = new BisulfiteDiploidSNPGenotypePriors();
		 
		 if(NAC.ptMode){
			 File fncw = new File(NAC.ocwf);
			 callableWindWriter = new bedObjectWriterImp(fncw);
			 String bedHeadLineWind = "chr\tstart\tend\taveMethyWind\tgchNumWind\tgchDepthWind\tgchCTdepthWind\tgchDotWind\tsigValue\taveMethyWindLinker\twcgMethyWind\twcgNumWind\n";
			 callableWindWriter.addHeader(bedHeadLineWind);
		 }
		 if(NAC.vm){
			 gchListForVerboseMode = new ArrayList<Double[]>();
		 }
		 summary = new ContextCondition();
		 rodIt = getToolkit().getRodDataSources().get(0);
	 }
	
	@Override
	public NDRCallContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
		// TODO Auto-generated method stub
		String cytosinePattern = "GCH-2";
		String wcgPattern = "WCG-2";
		//String hcgPattern = "HCG-2";
		double methyStatus = 0.36; //H1: 0.36; imr90:0.45
		double wcgMethyStatus = 0.80; //H1: 0.36; imr90:0.45
		//double hcgMethyStatus = 0.80;
		BisSNPUtils it = new BisSNPUtils(NAC);
		AlignmentContext stratifiedContexts = it.getFilteredAndStratifiedContexts(NAC, ref, context);

		NDRCallContext ncc = new NDRCallContext(stratifiedContexts, ref.getLocus()); // this make some biad for seq poor region's nuc window 
		
		if(stratifiedContexts == null){
			ncc.setCytosinePatternFlag(false);
			return ncc;
		}
			
		if(it.checkCytosineStatus(cytosinePattern, stratifiedContexts.getBasePileup(), tracker, ref, (BisulfiteDiploidSNPGenotypePriors) genotypePriors, NAC, methyStatus)){
			ncc.setCytosinePatternFlag(true);
		}
		//else if(it.checkCytosineStatus(hcgPattern, stratifiedContexts.getBasePileup(), tracker, ref, (BisulfiteDiploidSNPGenotypePriors) genotypePriors, NAC, hcgMethyStatus)){
		//	ncc.setHcgPatternFlag(true);
		else if(it.checkCytosineStatus(wcgPattern, stratifiedContexts.getBasePileup(), tracker, ref, (BisulfiteDiploidSNPGenotypePriors) genotypePriors, NAC, wcgMethyStatus)){
				ncc.setWcgPatternFlag(true);
		}
		//}
		
		
		
			//make fake context in chrM, so we know it immediately it is not Gch
			
			return ncc;
		//return NDRD_engine.calculateNDRscore(tracker, ref, context);
			
	}

	@Override
	public windowsObject reduce(NDRCallContext value, windowsObject windows) {
		// TODO Auto-generated method stub
		
		
		if ( value == null )
            return windows;

		if(windows.windowsPre.size() <= NAC.nucLinkerWindow - 1){
			if(windows.windowsPre.peekLast() != null){
				if(value.getLoc().discontinuousP(windows.windowsPre.peekLast().getLoc())){
					clearStatus(true);
				}
			}
			windows.windowsPre.offerLast(value);
	//		System.err.println(windows.windowsPre.size() + "\t" + windows.windowsMid.size() + "\t" + windows.windowsPost.size());
			return windows;
		}	
		else if(windows.windowsMid.size() <= NAC.nucPosWindow - 2){
			if(windows.windowsMid.peekLast() != null){
				if(value.getLoc().discontinuousP(windows.windowsMid.peekLast().getLoc())){
					clearStatus(true);
					
					windows.windowsPre.offerLast(value);
					
					return windows;
				}
			}
			windows.windowsMid.offerLast(value);
		//	System.err.println(windows.windowsPre.size() + "\t" + windows.windowsMid.size() + "\t" + windows.windowsPost.size());
			return windows;
		}
		else if(windows.windowsMid.size() == NAC.nucPosWindow - 1){
			if(windows.windowsMid.peekLast() != null){
				if(value.getLoc().discontinuousP(windows.windowsMid.peekLast().getLoc())){ // there is potential bias here for not enough reads in postWindow..
					clearStatus(true);
					
					windows.windowsPre.offerLast(value);
					
					return windows;
				}
			}
			windows.windowsMid.offerLast(value);
		//	System.err.println(windows.windowsPre.size() + "\t" + windows.windowsMid.size() + "\t" + windows.windowsPost.size());
		}
		else if(windows.windowsPost.size() <= NAC.nucLinkerWindow - 2){
			if(windows.windowsPost.peekLast() != null){
				if(value.getLoc().discontinuousP(windows.windowsPost.peekLast().getLoc())){
					
					
					winForceEndflag = true;
					
							if(checkNearRodLocation(windows.windowsMid.peekFirst().getLoc(),up,down, true, false) && checkNearRodLocation(windows.windowsMid.peekLast().getLoc(),up, down, true, false)){
								getNDRFromWindowsBySigTest(windows);
								addNDRtoWriter();
							}
							
						
					
					clearStatus(true);
					windows.windowsPre.offerLast(value);
					
					return windows;
				}
			}
			windows.windowsPost.offerLast(value);
		//	System.err.println(windows.windowsPre.size() + "\t" + windows.windowsMid.size() + "\t" + windows.windowsPost.size());
			return windows;
		}
		else if(windows.windowsPost.size() == NAC.nucLinkerWindow - 1){
			if(windows.windowsPost.peekLast() != null){
				if(value.getLoc().discontinuousP(windows.windowsPost.peekLast().getLoc())){
					
					
					winForceEndflag = true;
					
							if(checkNearRodLocation(windows.windowsMid.peekFirst().getLoc(),up,down, true, false) && checkNearRodLocation(windows.windowsMid.peekLast().getLoc(),up, down, true, false)){
								getNDRFromWindowsBySigTest(windows);
								addNDRtoWriter();
							}
						

					clearStatus(true);
					windows.windowsPre.offerLast(value);
					
					return windows;
				}
			}
			windows.windowsPost.offerLast(value);
		//	System.err.println(windows.windowsPre.size() + "\t" + windows.windowsMid.size() + "\t" + windows.windowsPost.size());
		}
		else{
			if(windows.windowsPost.peekLast() != null){
			//	if(windows.windowsMid.peekLast().getLoc().discontinuousP(windows.windowsPost.peekFirst().getLoc())){
				if(value.getLoc().discontinuousP(windows.windowsPost.peekLast().getLoc())){
					winForceEndflag = true;
					
							if(checkNearRodLocation(windows.windowsMid.peekFirst().getLoc(),up,down, true, false) && checkNearRodLocation(windows.windowsMid.peekLast().getLoc(),up, down, true, false)){
								getNDRFromWindowsBySigTest(windows);
								addNDRtoWriter();
							}
						

					clearStatus(true);
					windows.windowsPre.offerLast(value);
					
					return windows;
				}
			}
			windows.windowsPre.pollFirst();
			windows.windowsPre.offerLast(windows.windowsMid.pollFirst());
			windows.windowsMid.offerLast(windows.windowsPost.pollFirst());
			windows.windowsPost.offerLast(value);
			
		//	System.err.println(windows.windowsPre.size() + "\t" + windows.windowsMid.size() + "\t" + windows.windowsPost.size());
		}
		//summary.nWindowsVisited++;
		//if(NAC.largeLinker && !value.getCytosinePatternFlag() && !windows.windowsMid.peekLast().getCytosinePatternFlag() && !windows.windowsPre.peekLast().getCytosinePatternFlag()){
		//	return windows;
		//}
		//if(windows.windowsMid.peekFirst().getLoc().getStart()>= 16665 && windows.windowsMid.peekFirst().getLoc().getStart()<=16765)
		//	System.err.println(windows.windowsMid.peekFirst().getLoc().getStart() + "\t" + windows.windowsMid.peekLast().getLoc() + "\t" + winStartflag + "\t" + winEndflag + "\t" + winForceEndflag + "\trod: ");
		if(checkNearRodLocation(windows.windowsMid.peekFirst().getLoc(),up,down, true, false) && checkNearRodLocation(windows.windowsMid.peekLast().getLoc(),up, down, true, false)){
			//System.err.println(windows.windowsMid.peekFirst().getLoc().getStart() + "\t" + windows.windowsMid.peekLast().getLoc() + "\t" + winStartflag + "\t" + winEndflag + "\t" + winForceEndflag + "\trod: ");
			if(!outputFlag){
				
				//if(checkNearRodLocation(windows.windowsMid.peekFirst().getLoc(),1, 0, true)){
					winForceEndflag = true;
					//System.err.println(windows.windowsMid.peekFirst().getLoc().getStart());
				//}
				
				GenomeLoc tmp = windows.windowsMid.peekLast().getLoc();
				boolean rodNextExist = checkNearRodLocation(tmp,1,(NAC.nucPosWindow+2*NAC.nucLinkerWindow),false, false);
				getNDRFromWindowsBySigTest(windows);
			//	System.err.println(windows.windowsMid.peekFirst().getLoc().getStart() + "\t" + winStartflag + "\t" + winEndflag + "\t" + winForceEndflag + "\trod: " + rodNextExist);
				
			//	System.err.println(tmpGchNumWind + "\t" + tmpGchDepthWind + "\t" + tmpGchCTDepthWind + "\t" + tmpGchDotWind + "\t" + tmpGchMethyWind + "\t" + sigValueMem);
				addNDRtoWriter(!rodNextExist);
				if(winForceEndflag){
					clearStatus(!rodNextExist);
				}
			}
	
		}
		else{
			//GenomeLoc tmp = windows.windowsMid.peekLast().getLoc();
			//boolean rodNextExist = checkNearRodLocation(tmp,-1,(2*NAC.nucLinkerWindow),false);
			//clearStatus(!rodNextExist);
			outputFlag = false;
		}
		
		

		
		return windows;
		
	
	}

	@Override
	public windowsObject reduceInit() {
		// TODO Auto-generated method stub
		windows = new windowsObject();
		return windows;
	}
	
	@Override
	public windowsObject treeReduce(windowsObject lhWinds,
			windowsObject rhWinds) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public void onTraversalDone(windowsObject sum) {
//		logger.info(String.format("Visited windows                                %d", summary.nWindowsVisited));
		logger.info(String.format("Callable windows                                %d", summary.nWindowsCallable));
		
		//logger.info(String.format("Percentage of callable windows of all windows                                %.2f", summary.percentCallableWindowOfAll()));
		
		logger.info(String.format("Average GCH methylation in callable windows                                %.2f", summary.percentGchMethyOfCallableWindows()));
		//logger.info(String.format("Average HCG methylation in callable windows                                %.2f", summary.percentHcgMethyOfCallableWindows()));
		logger.info(String.format("Average WCG methylation in callable windows                                %.2f", summary.percentWcgMethyOfCallableWindows()));
		logger.info(String.format("Average GCH reads depth in callable windows                                %.2f", summary.percentGchDepthOfCallableWindows()));
		logger.info(String.format("Average GCH CT reads depth in callable windows                                %.2f", summary.percentGchCTDepthOfCallableWindows()));
		
		logger.info(String.format("Average GCH number in callable windows                                %.2f", summary.percentGchNumOfCallableWindows()));
		//logger.info(String.format("Average HCG number in callable windows                                %.2f", summary.percentHcgNumOfCallableWindows()));
		logger.info(String.format("Average WCG number in callable windows                                %.2f", summary.percentWcgNumOfCallableWindows()));
		logger.info(String.format("Average GCH data point in callable windows                                %.2f", summary.percentGchDotOfCallableWindows()));
		
	
		
			callableWindWriter.close();
		
	}
	
	public boolean checkNearRodLocation(GenomeLoc locus, int upstream, int downstream, boolean oneRod, boolean noDirection){
		GenomeLoc searchLoc = getToolkit().getGenomeLocParser().createGenomeLoc(locus.getContig(), locus.getStart() - Math.max(upstream, downstream)-1, locus.getStart() + Math.max(upstream, downstream));
		

		LocationAwareSeekableRODIterator locRodIt = rodIt.seek(searchLoc);
		//if(searchLoc.getStart()>=1027383 && searchLoc.getStart()<=1027483){
		//	SimpleBEDFeature bedtmp = (SimpleBEDFeature)locRodIt.seekForward(searchLoc).get(0).getUnderlyingObject();
		//	System.err.println("what is:" + locus.getStart() + "\t" + bedtmp.getStart() + "\t" + (bedtmp.getStart()-upstream) + "\t" + (bedtmp.getStart() + downstream) + "\t" + bedtmp.getStrand().toString());
		//}
    	// boolean locPre = locRodItPre.hasNext();
    	// boolean locPost = locRodItPost.hasNext();
		boolean rodExist = false;
    	 if(locRodIt.hasNext()){
    		RODRecordList rodList = locRodIt.seekForward(searchLoc);
    		
    		
    		if(rodList.size()<=1){
    			if(oneRod){
    				SimpleBEDFeature bed = (SimpleBEDFeature)locRodIt.seekForward(searchLoc).get(0).getUnderlyingObject();
    				
    				
    				//System.err.println(locus.getStart() + "\t" +bed.getStart() + "\t" + (downstream) + "\t" + (upstream) );
           		 	if((bed.getStrand()==Strand.POSITIVE && (locus.getStart() >= bed.getStart()-upstream && locus.getStart() < bed.getStart() + downstream)) || ((bed.getStrand()==Strand.NEGATIVE && (locus.getStart() > bed.getStart() - downstream && locus.getStart()<= bed.getStart() + upstream))) || (noDirection && (Math.abs(bed.getStart()-locus.getStart()) <= downstream || Math.abs(locus.getStart() - bed.getStart()) <= upstream))){
           		// 	System.err.println(locus.getStart() + "\t" + bed.getStart() + "\t" + (bed.getStart()-upstream) + "\t" + (bed.getStart() + downstream) + "\t" + bed.getStrand().toString());
           		 		rodExist=true; 
           		 	}
    			}
    			
    			
    		}
    		else{
    			if(oneRod){
    				SimpleBEDFeature bed = (SimpleBEDFeature)locRodIt.seekForward(searchLoc).get(0).getUnderlyingObject();
    				//System.err.println(locus.getStart() + "\t" +bed.getStart() + "\t" + (downstream) + "\t" + ( upstream) );
    				if((bed.getStrand()==Strand.POSITIVE && (locus.getStart() >= bed.getStart()-upstream && locus.getStart() < bed.getStart() + downstream)) || ((bed.getStrand()==Strand.NEGATIVE && (locus.getStart() > bed.getStart() - downstream && locus.getStart()<= bed.getStart() + upstream))) || (noDirection && (Math.abs(bed.getStart()-locus.getStart()) <= downstream || Math.abs(locus.getStart() - bed.getStart()) <= upstream))){
           		 		rodExist=true; 
           		 	}
    			}
    			else{
    				Iterator<GATKFeature> iter = locRodIt.seekForward(searchLoc).iterator();
    				while(iter.hasNext()){
    					SimpleBEDFeature bed = (SimpleBEDFeature)iter.next().getUnderlyingObject();
    					if((bed.getStrand()==Strand.POSITIVE && (locus.getStart() >= bed.getStart()-upstream && locus.getStart() < bed.getStart() + downstream)) || ((bed.getStrand()==Strand.NEGATIVE && (locus.getStart() > bed.getStart() - downstream && locus.getStart()<= bed.getStart() + upstream))) || (noDirection && (Math.abs(bed.getStart()-locus.getStart()) <= downstream || Math.abs(locus.getStart() - bed.getStart()) <= upstream))){
               		 		rodExist=true;
               		 		break;
               		 	    
               		 	}
    				}
    				
    			}
    			
       		 //SimpleBEDFeature bedPre = (SimpleBEDFeature)locRodItPre.seekForward(searchLocPre).get(0).getUnderlyingObject();
       		 //SimpleBEDFeature bedPost = (SimpleBEDFeature)locRodItPost.seekForward(searchLocPost).get(0).getUnderlyingObject();
    			//System.err.println(locus.getStart() + "\t" +bed.getStart());
    		}
    		 
    	 }

    	 rodIt.close(locRodIt);
    	// rodIt.close(locRodItPost);
		return rodExist;
	}
	
	
	
	// by statistics test
	public void getNDRFromWindowsBySigTest(windowsObject windows){
		
		
		windowsReturnObject objMid = getGchListFromWindow(windows.windowsMid);	
		windowsReturnObject objPre = getGchListFromWindow(windows.windowsPre);		
		windowsReturnObject objPost = getGchListFromWindow(windows.windowsPost);
		
		//System.err.println(windows.windowsMid.peekFirst().getLoc().getStart() + "\t" + windows.windowsMid.peekLast().getLoc().getStart());
		
		//System.err.println(objMid.numGchDot + "\t" + objMid.numValidGch + "\t" + objPre.numGchDot + "\t" + objPre.numValidGch + "\t" + objPost.numGchDot + "\t" + objPost.numValidGch);
		//something wrong here, callable window are not so many. 
		
		//if(objMid.numValidGch >= NAC.minGchNum){
		if(objMid.numGchDot >= NAC.minGchDotWindow && objMid.numValidGch >= NAC.minGchNum){	
			
			
			cytosineStat num = validateGch(windows.windowsMid.peekLast());
			int numC = 0;
			int numT = 0;
			int numO = 0;
			if(num != null){
				numC = num.numC;
				numT = num.numT;
				numO =num.numOther;
			}
			winChr = windows.windowsMid.peekFirst().getLoc().getContig();
			//look at windows callable or not
			if(!winEndflag){
				if(winStartflag && winInflag){
					if(num != null){
						if(num.type==cytosineType.GCH){
							tmpGchNumWind++;
							tmpGchCTDepthWind = tmpGchCTDepthWind + numC + numT;
							tmpGchDotWind = tmpGchDotWind + numC + numT;
							tmpGchNumCDotWind = tmpGchNumCDotWind + numC;
							tmpGchDepthWind = tmpGchDepthWind + numC + numT + numO;
							tmpGchMethyWind += (double)numC/(double)(numC + numT);
						}
						else if(num.type==cytosineType.WCG){
							tmpWcgNumWind++;
							tmpWcgMethyWind += (double)numC/(double)(numC + numT);
						}
						
					}
					if(winForceEndflag){
						winEndflag = true;
						winEndCor = windows.windowsMid.peekLast().getLoc().getStart();
					}
				}
				else{
					//if(objPre.numValidGch >= NAC.minGchNumLinkerWindow){
					if(objPre.numGchDot >= NAC.minGchDotLinkerWindow && objPre.numValidGch >= NAC.minGchNumLinkerWindow){
						winPreLinkerflag = true;
						winStartflag = true;
						winInflag = true;
						winStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
						tmpGchNumWind = objMid.numValidGch;
						tmpGchCTDepthWind = objMid.sumGchCTDepth;
						tmpGchDotWind = objMid.numGchDot;
						tmpGchNumCDotWind = objMid.numCOfGchDot;
						tmpGchDepthWind = objMid.sumGchDepth;
						tmpGchMethyWind = objMid.sumGchMethy;
						tmpWcgMethyWind = objMid.sumWcgMethy;
						tmpWcgNumWind = objMid.numValidWcg;
						double sigValue;
						if(NAC.test.equalsIgnoreCase("binomialTest")){
							sigValue = getBinomialSigTest(objMid.numCOfGchDot, objMid.numGchDot, objPre.sumGchMethy/objPre.numValidGch);
						}
						else{
							sigValue = getSigTest(objMid.dot, objPre.dot);
						}
							
						
						
							if(this.sigValueMem >= sigValue || this.sigValueMem == -1){
								this.sigValueMem = sigValue;
								this.tmpGchMethyWindLinker = objPre.sumGchMethy/objPre.numValidGch;
							}
								
						
						if(winForceEndflag){
							winEndflag = true;
							winEndCor = windows.windowsMid.peekLast().getLoc().getStart();
						}
					}
					//if(objPost.numValidGch >= NAC.minGchNumLinkerWindow){
					if(objPost.numGchDot >= NAC.minGchDotLinkerWindow && objPost.numValidGch >= NAC.minGchNumLinkerWindow){
						winPostLinkerflag = true;
						
						if(!winPreLinkerflag){
							winStartflag = true;
							winInflag = true;
							winStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
							tmpGchNumWind = objMid.numValidGch;
							tmpGchCTDepthWind = objMid.sumGchCTDepth;
							tmpGchDotWind = objMid.numGchDot;
							tmpGchNumCDotWind = objMid.numCOfGchDot;
							tmpGchDepthWind = objMid.sumGchDepth;
							tmpGchMethyWind = objMid.sumGchMethy;
							tmpWcgMethyWind = objMid.sumWcgMethy;
							tmpWcgNumWind = objMid.numValidWcg;
							double sigValue;
							if(NAC.test.equalsIgnoreCase("binomialTest")){
								sigValue = getBinomialSigTest(objMid.numCOfGchDot, objMid.numGchDot, objPost.sumGchMethy/objPost.numValidGch);
							}
							else{
								sigValue = getSigTest(objMid.dot, objPost.dot);
							}
							if(winInflag){
								if(this.sigValueMem >= sigValue || this.sigValueMem == -1){
									this.sigValueMem = sigValue;
									this.tmpGchMethyWindLinker = objPost.sumGchMethy/objPost.numValidGch;
								}
									
							}
							if(winForceEndflag){
								winEndflag = true;
								winEndCor = windows.windowsMid.peekLast().getLoc().getStart();
							}
						}
						
						
					}
					
				}
			}
			
			
		}
		else{
			if(winStartflag && winInflag){
				winEndflag = true;
			//	winStartflag = false;
			//	winInflag = false;
				winEndCor = windows.windowsMid.peekLast().getLoc().getStart();

			}
			
			
		}
		
		
	}

	
	public class windowsObject {
		public LinkedList<NDRCallContext> windowsMid = null;
		public LinkedList<NDRCallContext> windowsPre = null;
		public LinkedList<NDRCallContext> windowsPost = null;
		public boolean windowsPreContinuous = true;
		public boolean windowsPostContinuous = true;
		public windowsObject(){
			windowsMid = new LinkedList<NDRCallContext>();
			windowsPre = new LinkedList<NDRCallContext>();
			windowsPost = new LinkedList<NDRCallContext>();
		}
	
	}
	
	public windowsReturnObject getGchListFromWindow(LinkedList<NDRCallContext> window){
		double sumGchMethy = 0;
		int numValidGch = 0;
		int numGchDot = 0;
		int numCOfGchDot = 0;
		int sumGchCTDepth = 0;
		int sumGchDepth = 0;
		double sumWcgMethy = 0;
		int numValidWcg = 0;
		//ArrayList<Double> dot = new ArrayList<Double>();
		ArrayList<Double> dot = new ArrayList<Double>();
		Iterator<NDRCallContext> itContext = window.iterator();
		while(itContext.hasNext()){
			NDRCallContext tmpContext = itContext.next();
			cytosineStat num = validateGch(tmpContext);
			
			if(num != null){
				int numC = num.numC;
				int numT = num.numT;
				int numO = num.numOther;
				if(num.type==cytosineType.GCH){
					numValidGch++;
					sumGchCTDepth += (numC + numT);
					numGchDot += (numC + numT);
					numCOfGchDot += numC;
					sumGchDepth += (numC + numT + numO);
					sumGchMethy += (double)numC/(double)(numC + numT);
					dot.add((double)numC/(double)(numC + numT));
				}
				else if(num.type==cytosineType.WCG){
					sumWcgMethy += (double)numC/(double)(numC + numT);
					numValidWcg++;
				}
				
				//for(int i=numC; i>=0;i--){
				//	dot.add(1);
				//}
				//for(int i=numT; i>=0;i--){
				//	dot.add(0);
				//}
				
			 }	
			//	System.out.println("loc: " + tmpContext.getLoc().getStart() + "\tGchMethy: " + (double)numC/(double)(numC + numT) + "\tnumC: " + numC + "\tnumT: " + numT);
			
		}
		windowsReturnObject obj = new windowsReturnObject(numValidGch, sumGchDepth, sumGchCTDepth, sumGchMethy, numGchDot,numCOfGchDot, dot, sumWcgMethy, numValidWcg);
		return obj;
	}

	public class windowsReturnObject {
		public int numValidGch;
		public int numGchDot;
		public int numCOfGchDot;
		public int sumGchCTDepth;
		public int sumGchDepth;
		public double sumGchMethy;
		public int numValidWcg;
		public double sumWcgMethy;
		public ArrayList<Double> dot;
		public windowsReturnObject(int numValidGch,int sumGchDepth, int sumGchCTDepth,double sumGchMethy, int numGchDot,int numCOfGchDot, ArrayList<Double> dot, double sumWcgMethy,int numValidWcg){
			this.numValidGch = numValidGch;
			this.sumGchCTDepth = sumGchCTDepth;
			this.sumGchDepth = sumGchDepth;
			this.sumGchMethy = sumGchMethy;
			this.numGchDot = numGchDot;
			this.numCOfGchDot = numCOfGchDot;
			this.dot = dot;
			this.numValidWcg = numValidWcg;
			this.sumWcgMethy = sumWcgMethy;
		}
	}
	
	public cytosineStat validateGch(NDRCallContext ncc){
		cytosineStat returnObj = new cytosineStat();
		//int[] num = new int[3]; //0: number of C; 1: number of T; 2: number of other bases;
		//num[0] = 0;
		//num[1] = 0;
		//num[2] = 0;
		if(ncc.getCytosinePatternFlag()){
			returnObj.type = cytosineType.GCH;
		}
		//else if(ncc.getHcgPatternFlag()){
		//	returnObj.type = cytosineType.HCG;
		else if(ncc.getWcgPatternFlag()){
				returnObj.type = cytosineType.WCG;
			}
		//}
		else{
			return null;
		}
			
		int numC = 0;
		int numT = 0;
		int numO = 0;
		for( PileupElement p : ncc.getRealContext().getBasePileup()){
			SAMRecord samRecord = p.getRead();
			if(samRecord.getDuplicateReadFlag()){ //get rid of duplicate reads
            	continue;
            }
        	int offset = p.getOffset();
        	if(offset < 0)//is deletion
        		continue;
        	boolean paired = samRecord.getReadPairedFlag();
        	if(paired){
        		try {
					samRecord = (SAMRecord) p.getRead().clone();
				} catch (CloneNotSupportedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
	        	boolean secondOfPair = samRecord.getSecondOfPairFlag();

	        	if (samRecord.getNotPrimaryAlignmentFlag())
				{
					continue;
				}
				
				// Inverted dups, count only one end
				if (samRecord.getAlignmentStart() == samRecord.getMateAlignmentStart() && samRecord.getReadNegativeStrandFlag() == samRecord.getMateNegativeStrandFlag())
				{
					if (samRecord.getSecondOfPairFlag()) continue;
   				}
	        	if (paired  && !NAC.USE_BADLY_MATED_READS && !samRecord.getProperPairFlag())
				{
					continue;
				}
	        	
	        	
	        	if(secondOfPair){	        		
		        	samRecord.setReadNegativeStrandFlag(!samRecord.getReadNegativeStrandFlag());        		
	        	}
	        	
        	}
			
        	boolean negStrand = samRecord.getReadNegativeStrandFlag();
		//	int alignmentS = samRecord.getAlignmentStart();
		//	int	onRefCoord = (negStrand) ? samRecord.getUnclippedEnd() : alignmentS;
			
			
			if(((GATKSAMRecord)p.getRead()).isGoodBase(p.getOffset())){
				if(negStrand){
					if(p.getBase()==BaseUtils.G){
						numC++;
					}
					else if(p.getBase()==BaseUtils.A){
						numT++;
					}
					else{
						numO++;
					}
					
				}
				else{
					if(p.getBase()==BaseUtils.C){
						numC++;
					}
					else if(p.getBase()==BaseUtils.T){
						numT++;
					}
					else{
						numO++;
					}
				}
			}
			
		}
		if((numC + numT) >= NAC.minCTDepth && (numC + numT + numO) >= NAC.minDepth){
			returnObj.numC = numC;
			returnObj.numT = numT;
			returnObj.numOther = numO;
			return returnObj;
		}
		else{
			return null;
		}
		
	}
	
	//using MannWhitneyTest/Wilcox rank sum test
	public double getSigTest(ArrayList<Double> dotMid, ArrayList<Double> dotAdj){
		
		double[] mid = new double[dotMid.size()];
		double[] adj = new double[dotAdj.size()];
		for(int i = 0; i < dotMid.size();i++){
			mid[i] = dotMid.get(i);
		}
		for(int i = 0; i < dotAdj.size();i++){
			adj[i] = dotAdj.get(i);
		}
		if(NAC.test.equalsIgnoreCase("rankSumTest")){
			MannWhitneyTest test = new MannWhitneyTest(mid, adj, H1.GREATER_THAN);
			return test.getSP();	
		}
		else if(NAC.test.equalsIgnoreCase("tTest")){
			TwoSampleTtest test = new TwoSampleTtest(mid, adj, H1.GREATER_THAN, NAC.samVar);
			return test.getSP();	
		}
		else if(NAC.test.equalsIgnoreCase("ksTest")){
			SmirnovTest test = new SmirnovTest(mid, adj, H1.GREATER_THAN);
			return test.getSP();	
		}
		else{
			System.err.println("Wrong test name!");
			System.exit(1);
		}
		return Double.POSITIVE_INFINITY;
		
		//
		
	}
	
	public double getBinomialSigTest(int k, int n, double pValue){
		
		//double p;
		//if(k==0){
		//	p = Probability.binomial(k,n,pValue);
		//}
		//else{
		//	p = Probability.binomial(k,n,pValue)-Probability.binomial(k-1,n,pValue);
		//}
		//System.err.println(k + "\t" + n + "\t" + pValue + "\t" + p);
		BinomialDistributionImpl binomial = new BinomialDistributionImpl(n, pValue);
		
		return binomial.probability(k);
		//return p;
		//
		
	}
	
	
	public void addNDRtoWriter(){
		addNDRtoWriter(true);
	}
	
	public void addNDRtoWriter(boolean clearWind){
		//boolean outputFlag = false;
		//if(NDREndflag && tmpGchNumWindNDR  >= NAC.minGchNum && tmpGchCTDepthWindNDR/tmpGchNumWindNDR >= NAC.minCTDepth && tmpGchNumWindLinker  >= NAC.minGchNumLinkerWindow && tmpGchCTDepthWindLinker/tmpGchNumWindLinker >= NAC.minCTDepthLinkerWindow){
		

		
			//if(winEndflag && tmpGchNumWind  >= (NAC.minGchNum) && tmpGchCTDepthWind/tmpGchNumWind >= NAC.minCTDepth && this.sigValueMem != -1){
			if(winEndflag && winStartflag && winInflag){
				ArrayList<Object> valuesWind = new ArrayList<Object>();
				
				valuesWind.add(tmpGchMethyWind/tmpGchNumWind);
				valuesWind.add(tmpGchNumWind);
				valuesWind.add(tmpGchDepthWind/tmpGchNumWind);
				valuesWind.add(tmpGchCTDepthWind/tmpGchNumWind);
				valuesWind.add(tmpGchDotWind);	
				valuesWind.add(this.sigValueMem);
				//valuesWind.add(tmpGchNumCDotWind);
				valuesWind.add(tmpGchMethyWindLinker);
				valuesWind.add(tmpWcgMethyWind/tmpWcgNumWind);
				valuesWind.add(tmpWcgNumWind);
				
				
				bedObject bedLineWind = new bedObject(winChr, winStartCor, winEndCor, valuesWind); //chr, start, end, aveMethyWind, gchNumWind, gchCTdepthWind
				callableWindWriter.add(bedLineWind);
				summary.nWindowsCallable++;
				summary.sumGchDepthWind += tmpGchDepthWind;
				summary.sumGchCTDepthWind += tmpGchCTDepthWind;
				summary.sumGchMethyWindowsCalledConfidently += tmpGchMethyWind;
				summary.sumGchNumInWindCalledConfidently += tmpGchNumWind;
				summary.sumGchDotWind += tmpGchDotWind;
				
				summary.sumWcgNumInWindCalledConfidently += tmpWcgNumWind;
				summary.sumWcgMethyWindowsCalledConfidently += tmpWcgMethyWind;
				
				outputFlag=true;
				clearStatus(clearWind);
			}
			

		

	}
	
	private void clearStatus(boolean clearWind){
		if(clearWind){
			windows.windowsPre.clear();
			windows.windowsMid.clear();
			windows.windowsPost.clear();
		}
		
		clearFlagStatus();

	}
	
	
	private void clearFlagStatus(){
		
		winStartflag = false;	
		winEndflag = false;	
		winForceEndflag = false;
		winInflag = false;	
		winStartCor = -1;			
		winEndCor = -1;
		winChr = null;
		winPreLinkerflag = false;
		winPostLinkerflag = false;
		sigValueMem = -1;
		
		tmpGchCTDepthWind = 0;
		tmpGchDepthWind = 0;
		tmpGchNumWind = 0;
		tmpGchMethyWind = 0;
		tmpGchDotWind = 0;
		tmpGchNumCDotWind = 0;
		
		tmpWcgNumWind = 0;
		tmpWcgMethyWind = 0;
		
		
	}
	
	
	private class cytosineStat{
		int numC;
		int numT;
		int numOther;
		cytosineType type;
	}
	private enum cytosineType{
		GCH,
		HCG,
		WCG
	}
	
	
	private boolean passWindowFilter(NDRCallContext value){
		if(NAC.ptMode){
			if ( value == null ){
				return false;
			}
			else{
				if(value.hasRealContext()){
					if(value.getRealContext().getBasePileup().size() > NAC.minDepth){
						if(value.getCytosinePatternFlag()){
							if(validateGch(value) == null){
								return false;
							}
						
						}
						
					}
					else{
						return false;
					}
				}
				else{
					return false;
				}
			}
		}
		return true;
	}



}
