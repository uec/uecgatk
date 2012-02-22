package edu.usc.epigenome.uecgatk.YapingWalker;

/*
 * program measure NDR from the lef first position with sig value to the right first position with sig value. so NDR position would have a little shift (may be 15bp) to the left position. 
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;


import jsc.independentsamples.MannWhitneyTest;
import jsc.independentsamples.SmirnovTest;
import jsc.independentsamples.TwoSampleTtest;
import jsc.tests.H1;

import net.sf.samtools.SAMRecord;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentReadFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

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
import edu.usc.epigenome.uecgatk.YapingWalker.NDRdetectWalker.windowsObject;
import edu.usc.epigenome.uecgatk.YapingWalker.NDRargumentCollection;

@BAQMode(QualityMode = BAQ.QualityMode.OVERWRITE_QUALS, ApplicationTime = BAQ.ApplicationTime.ON_INPUT)
@ReadFilters( {UnmappedReadFilter.class, BadMateFilter.class, NotPrimaryAlignmentReadFilter.class, DuplicateReadFilter.class} )
@Reference(window=@Window(start=-300,stop=300))
@By(DataSource.REFERENCE)
@Downsample(by=DownsampleType.NONE)
public class NDRdetectWalker extends LocusWalker<NDRCallContext,windowsObject> implements
		TreeReducible<windowsObject> {

	@ArgumentCollection private static NDRargumentCollection NAC = new NDRargumentCollection();
	
	
    protected bedObjectWriterImp writer = null;
    
    protected bedObjectWriterImp callableWindWriter = null;
	
	//private NDRdetectionEngine NDRD_engine = null;
	
	private windowsObject windows = null;
	
	private GenotypePriors genotypePriors;
	
	private ContextCondition summary = null;
	
	private boolean NDRStartflag = false;
	
	private boolean NDREndflag = false;
	
	private boolean NDRInflag = false;
	
	private int NDRStartCor = -1;
			
	private int NDREndCor = -1;
	
	private boolean NDRPreLinkerflag = false;
	
	private boolean winPreLinkerflag = false;
	
	private boolean NDRPostLinkerflag = false;
	
	private boolean winPostLinkerflag = false;
	
	private boolean winStartflag = false;
	
	private boolean winEndflag = false;
	
	private boolean winInflag = false;
	
	private boolean winForceEndflag = false;
	
	private boolean winPassFilterflag = true;
	
	private String winChr = null;
	
	private int winStartCor = -1;
			
	private int winEndCor = -1;
	
	private int tmpGchCTDepthWind = 0;
	
	private int tmpGchNumWind = 0;
	
	private int tmpGchDepthWind = 0;
	
	private int tmpGchDotWind = 0;
	

	
	private double tmpGchMethyWind = 0;
	
	private int tmpGchCTDepthWindLinker = 0;
	
	private int tmpGchDepthWindLinker = 0;
	
	private int tmpGchNumWindLinker = 0;
	
	private int tmpGchDotWindLinker = 0;
	
	private double tmpGchMethyWindLinker = 0;
	
	private int tmpGchNumWindNDR = 0;
	
	private int tmpGchCTDepthWindNDR = 0;
	
	private int tmpGchDepthWindNDR = 0;
	
	private int tmpGchDotWindNDR = 0;
	
	private double tmpGchMethyWindNDR = 0;
	
	private double sigValueMem = -1;
	
	private ArrayList<Double[]> gchListForVerboseMode = null;
	
	public NDRdetectWalker() {
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
	        double percentGchCTDepthOfCallableWindows() { return (double)sumGchCTDepthWind/(nWindowsCallable * percentGchNumOfCallableWindows());}
	        double percentGchDepthOfCallableWindows() { return (double)sumGchDepthWind/(nWindowsCallable * percentGchNumOfCallableWindows());}
	        double percentGchDotOfCallableWindows() { return (double)sumGchDotWind/(nWindowsCallable);}
	        double percentGchNumOfCallableWindows() { return (double)sumGchNumInWindCalledConfidently/(nWindowsCallable);}
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
		 File fn = new File(NAC.outFile);
		 writer = new bedObjectWriterImp(fn);
		 String bedHeadLine = "chr\tstart\tend\taveMethyNDR\tgchNumNDR\tgchDepthNDR\tgchCTdepthNDR\tgchDotNDR\taveMethyLinker\tgchNumLinker\tgchDepthLinker\tgchCTdepthLinker\tgchDotLinker\tsigValue\n";
		 writer.addHeader(bedHeadLine);
		 if(NAC.ptMode){
			 File fncw = new File(NAC.ocwf);
			 callableWindWriter = new bedObjectWriterImp(fncw);
			 String bedHeadLineWind = "chr\tstart\tend\taveMethyWind\tgchNumWind\tgchDepthWind\tgchCTdepthWind\tgchDotWind\tsigValue\n";
			 callableWindWriter.addHeader(bedHeadLineWind);
		 }
		 if(NAC.vm){
			 gchListForVerboseMode = new ArrayList<Double[]>();
		 }
		 summary = new ContextCondition();
	 }
	
	@Override
	public NDRCallContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
		// TODO Auto-generated method stub
		String cytosinePattern = "GCH-2";
		double methyStatus = 0.36; //H1: 0.36; imr90:0.45
		BisSNPUtils it = new BisSNPUtils(NAC);
		AlignmentContext stratifiedContexts = it.getFilteredAndStratifiedContexts(NAC, ref, context);
		
		NDRCallContext ncc = new NDRCallContext(stratifiedContexts, ref.getLocus()); // this make some biad for seq poor region's nuc window 
		
		if(stratifiedContexts == null){
			ncc.setCytosinePatternFlag(false);
			return ncc;
		}
			
		
		if(it.checkCytosineStatus(cytosinePattern, stratifiedContexts.getBasePileup(), tracker, ref, (BisulfiteDiploidSNPGenotypePriors) genotypePriors, NAC, methyStatus)){
			//System.out.println(ref.getLocus().getStart());
			 ncc.setCytosinePatternFlag(true);
			 return ncc;
		}	
		else{
			//make fake context in chrM, so we know it immediately it is not Gch
			ncc.setCytosinePatternFlag(false);
			return ncc;
		}
		//return NDRD_engine.calculateNDRscore(tracker, ref, context);
			
	}

	@Override
	public windowsObject reduce(NDRCallContext value, windowsObject windows) {
		// TODO Auto-generated method stub
		
		
		if ( value == null )
            return windows;

	//	System.err.println(windows.windowsPre.size() + "\t" + windows.windowsMid.size() + "\t" + windows.windowsPost.size());
		
		if(windows.windowsPre.size() <= NAC.nucLinkerWindow - 1){
			if(windows.windowsPre.peekLast() != null){
				if(value.getLoc().discontinuousP(windows.windowsPre.peekLast().getLoc())){
					clearStatus();
				}
			}
			windows.windowsPre.offerLast(value);
			return windows;
		}	
		else if(windows.windowsMid.size() <= NAC.nucPosWindow - 2){
			if(windows.windowsMid.peekLast() != null){
				if(value.getLoc().discontinuousP(windows.windowsMid.peekLast().getLoc())){
					clearStatus();
					winPassFilterflag = true;
					windows.windowsPre.offerLast(value);
					
					return windows;
				}
			}
			windows.windowsMid.offerLast(value);
			if(NAC.ptMode){
				if(!passWindowFilter(value)){
					//addNDRtoWriter();
					winPassFilterflag = false;
					//clearStatus();
					//return windows;
				}
			}
			return windows;
		}
		else if(windows.windowsMid.size() == NAC.nucPosWindow - 1){
			if(windows.windowsMid.peekLast() != null){
				if(value.getLoc().discontinuousP(windows.windowsMid.peekLast().getLoc())){ // there is potential bias here for not enough reads in postWindow..
					clearStatus();
					winPassFilterflag = true;
					windows.windowsPre.offerLast(value);
					
					return windows;
				}
			}
			windows.windowsMid.offerLast(value);
			if(NAC.ptMode){
				if(!passWindowFilter(value)){
					//addNDRtoWriter();
					winPassFilterflag = false;
					//clearStatus();
					//return windows;
				}
			}
		}
		else if(windows.windowsPost.size() <= NAC.nucLinkerWindow - 2){
			if(windows.windowsPost.peekLast() != null){
				if(value.getLoc().discontinuousP(windows.windowsPost.peekLast().getLoc())){
					if(NAC.ptMode && !winPassFilterflag){
						winPassFilterflag = true;
						clearStatus();
						windows.windowsPre.offerLast(value);
						return windows;
					}
					
					winForceEndflag = true;
					if(NAC.statTest){
						if(NAC.randomNoise){
							getNDRFromWindowsBySigTestWithRandomNoise(windows);
						}
						else{
							getNDRFromWindowsBySigTest(windows);
						}
					}
					addNDRtoWriter();
					clearStatus();
					windows.windowsPre.offerLast(value);
					
					return windows;
				}
			}
			windows.windowsPost.offerLast(value);
			
			return windows;
		}
		else if(windows.windowsPost.size() == NAC.nucLinkerWindow - 1){
			if(windows.windowsPost.peekLast() != null){
				if(value.getLoc().discontinuousP(windows.windowsPost.peekLast().getLoc())){
					if(NAC.ptMode && !winPassFilterflag){
						winPassFilterflag = true;
						clearStatus();
						windows.windowsPre.offerLast(value);
						return windows;
					}
					
					winForceEndflag = true;
					if(NAC.statTest){
						if(NAC.randomNoise){
							getNDRFromWindowsBySigTestWithRandomNoise(windows);
						}
						else{
							getNDRFromWindowsBySigTest(windows);
						}
					}
					addNDRtoWriter();
					clearStatus();
					windows.windowsPre.offerLast(value);
					
					return windows;
				}
			}
			windows.windowsPost.offerLast(value);
			
		}
		else{
			if(windows.windowsPost.peekLast() != null){
			//	if(windows.windowsMid.peekLast().getLoc().discontinuousP(windows.windowsPost.peekFirst().getLoc())){
				if(value.getLoc().discontinuousP(windows.windowsPost.peekLast().getLoc())){
					winForceEndflag = true;
					if(NAC.statTest){
						if(NAC.randomNoise){
							getNDRFromWindowsBySigTestWithRandomNoise(windows);
						}
						else{
							getNDRFromWindowsBySigTest(windows);
						}
					}
					addNDRtoWriter();
					clearStatus();
					windows.windowsPre.offerLast(value);
					
					return windows;
				}
			}
			windows.windowsPre.pollFirst();
			windows.windowsPre.offerLast(windows.windowsMid.pollFirst());
			windows.windowsMid.offerLast(windows.windowsPost.pollFirst());
			windows.windowsPost.offerLast(value);
			if(NAC.ptMode){
				if(!winPassFilterflag){
					
					winPassFilterflag = true;
					clearStatus();
					return windows;
				}
				if(!passWindowFilter(windows.windowsPost.peekFirst())){
					winForceEndflag = true;
					getNDRFromWindowsBySigTest(windows);
					addNDRtoWriter();
					clearStatus();
					return windows;
				}
			}
			
		}
		//summary.nWindowsVisited++;
		if(NAC.largeLinker && !value.getCytosinePatternFlag() && !windows.windowsMid.peekLast().getCytosinePatternFlag() && !windows.windowsPre.peekLast().getCytosinePatternFlag()){
			return windows;
		}
		
		if(NAC.statTest){
		//	if(NAC.randomNoise){
		//		getNDRFromWindowsBySigTestWithRandomNoise(windows);
		//	}
		//	else{
		//	System.err.println("ok");	
			getNDRFromWindowsBySigTest(windows);
		//	}
			
		}
		else{
			getValueFromWindow(windows.windowsMid);
		}
	//	
		//GenomeLoc locMidInWind = window.get(NAC.nucPosWindow/2).getLocation();
		addNDRtoWriter();
		
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
		logger.info(String.format("Confidantly called NDR windows                                %d", summary.nNDRWindowsCalledConfidently));
		//logger.info(String.format("Percentage of callable windows of all windows                                %.2f", summary.percentCallableWindowOfAll()));
		logger.info(String.format("Percentage of Confidantly called NDR windows of callable windows                                %.2f", summary.percentNDRWindowOfCallableWindows()));
		logger.info(String.format("Average GCH methylation in callable windows                                %.2f", summary.percentGchMethyOfCallableWindows()));
		logger.info(String.format("Average GCH reads depth in callable windows                                %.2f", summary.percentGchDepthOfCallableWindows()));
		logger.info(String.format("Average GCH CT reads depth in callable windows                                %.2f", summary.percentGchCTDepthOfCallableWindows()));
		
		logger.info(String.format("Average GCH number in callable windows                                %.2f", summary.percentGchNumOfCallableWindows()));
		logger.info(String.format("Average GCH data point in callable windows                                %.2f", summary.percentGchDotOfCallableWindows()));
		
		logger.info(String.format("Average GCH methylation in NDR windows                                %.2f", summary.percentGchMethyOfNDRWindowsCalledConfidently()));
		logger.info(String.format("Average GCH reads depth in NDR windows                                %.2f", summary.percentGchDepthOfNDRWindows()));
		logger.info(String.format("Average GCH CT reads depth in NDR windows                                %.2f", summary.percentGchCTDepthOfNDRWindows()));
		logger.info(String.format("Average GCH number in NDR windows                                %.2f", summary.percentGchNumOfNDRWindows()));
		logger.info(String.format("Average GCH data point in NDR windows                                %.2f", summary.percentGchDotOfNDRWindows()));
		
		logger.info(String.format("Average GCH methylation in NDR linker windows                                %.2f", summary.percentGchMethyOfCallableWindowsLinker()));
		logger.info(String.format("Average GCH reads depth in NDR linker windows                                %.2f", summary.percentGchDepthOfCallableWindowsLinker()));
		logger.info(String.format("Average GCH CT reads depth in NDR linker windows                                %.2f", summary.percentGchCTDepthOfCallableWindowsLinker()));
		logger.info(String.format("Average GCH number in NDR linker windows                                %.2f", summary.percentGchNumOfCallableWindowsLinker()));
		logger.info(String.format("Average GCH data point in NDR linker windows                                %.2f", summary.percentGchDotOfCallableWindowsLinker()));
		
		writer.close();
		if(NAC.ptMode){
			callableWindWriter.close();
		}
	}
	
	// by statistics test comparing with fake region simulate whole genome level GCH methylation
	public void getNDRFromWindowsBySigTestWithRandomNoise(windowsObject windows){
		
		
	}
	
	
	// by statistics test
	public void getNDRFromWindowsBySigTest(windowsObject windows){
		
		
		windowsReturnObject objMid = getGchListFromWindow(windows.windowsMid);	
		windowsReturnObject objPre = getGchListFromWindow(windows.windowsPre);		
		windowsReturnObject objPost = getGchListFromWindow(windows.windowsPost);
		
		
		
		
		//something wrong here, callable window are not so many. 
		
		//if(objMid.numValidGch >= NAC.minGchNum){
		if(objMid.numGchDot >= NAC.minGchDotWindow && objMid.numValidGch >= NAC.minGchNum){	
			
			
			int[] num = validateGch(windows.windowsMid.peekLast());
			int numC = 0;
			int numT = 0;
			int numO = 0;
			if(num != null){
				numC = num[0];
				numT = num[1];
				numO = num[2];
			}
			winChr = windows.windowsMid.peekFirst().getLoc().getContig();
			//look at windows callable or not
			if(!winEndflag){
				if(winStartflag && winInflag){
					if(num != null){
						tmpGchNumWind++;
						tmpGchCTDepthWind = tmpGchCTDepthWind + numC + numT;
						tmpGchDotWind = tmpGchDotWind + numC + numT;
						tmpGchDepthWind = tmpGchDepthWind + numC + numT + numO;
						tmpGchMethyWind += (double)numC/(double)(numC + numT);
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
						tmpGchDepthWind = objMid.sumGchDepth;
						tmpGchMethyWind = objMid.sumGchMethy;
						double sigValue = getSigTest(objMid.dot, objPre.dot);
						if(!NDRInflag && winInflag){
							if(this.sigValueMem >= sigValue || this.sigValueMem == -1){
								this.sigValueMem = sigValue;
							}
								
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
							tmpGchDepthWind = objMid.sumGchDepth;
							tmpGchMethyWind = objMid.sumGchMethy;
							double sigValue = getSigTest(objMid.dot, objPost.dot);
							if(!NDRInflag && winInflag){
								if(this.sigValueMem >= sigValue || this.sigValueMem == -1){
									this.sigValueMem = sigValue;
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
			
			//look at NDR
			if(!NDREndflag){
				if(NDRStartflag && NDRInflag){
					if(num != null){
						tmpGchNumWindNDR++;
						tmpGchCTDepthWindNDR = tmpGchCTDepthWindNDR + numC + numT;
						tmpGchDotWindNDR = tmpGchDotWindNDR + numC + numT;
						tmpGchDepthWindNDR = tmpGchDepthWindNDR + numC + numT + numO;
						tmpGchMethyWindNDR += (double)numC/(double)(numC + numT);
						if(NAC.vm){
							Double[] tmpGch = new Double[2];
							tmpGch[0] = (double) windows.windowsMid.peekLast().getLoc().getStart();
							tmpGch[1] =	(double)numC/(double)(numC + numT);
							gchListForVerboseMode.add(tmpGch);
						}
					}
					//if(objPost.numValidGch >= NAC.minGchNumLinkerWindow){
					if(objPost.numGchDot >= NAC.minGchDotLinkerWindow && objPost.numValidGch >= NAC.minGchNumLinkerWindow){
						double sigValue = getSigTest(objMid.dot, objPost.dot);
						if(sigValue < NAC.sigValue){
							NDREndflag = true;
							NDRPostLinkerflag = true;
							NDREndCor = windows.windowsMid.peekLast().getLoc().getStart();
							if(NDRPreLinkerflag){
								tmpGchNumWindLinker += objPost.numValidGch;
								tmpGchCTDepthWindLinker += objPost.sumGchCTDepth;
								tmpGchDotWindLinker += objPost.numGchDot;
								tmpGchDepthWindLinker += objPost.sumGchDepth;
								tmpGchMethyWindLinker += objPost.sumGchMethy;
							}
							else{
								
								sigValueMem = sigValue;

								tmpGchNumWindLinker = objPost.numValidGch;
								tmpGchCTDepthWindLinker = objPost.sumGchCTDepth;
								tmpGchDotWindLinker = objPost.numGchDot;
								tmpGchDepthWindLinker = objPost.sumGchDepth;
								tmpGchMethyWindLinker = objPost.sumGchMethy;
								
							}	
							winStartflag = true;
							winInflag = true;
							winEndflag = true;
							winEndCor = NDREndCor;
							winStartCor = NDRStartCor;
							tmpGchNumWind = tmpGchNumWindNDR;
							tmpGchCTDepthWind = tmpGchCTDepthWindNDR;
							tmpGchDotWind = tmpGchDotWindNDR;
							tmpGchDepthWind = tmpGchDepthWindNDR;
							tmpGchMethyWind = tmpGchMethyWindNDR;
						}
						
						
					}
					
					if(winForceEndflag){
						NDREndCor = windows.windowsMid.peekLast().getLoc().getStart();
						NDREndflag = true;
						winStartflag = true;
						winInflag = true;
						winEndflag = true;
						winEndCor = NDREndCor;
						winStartCor = NDRStartCor;
						tmpGchNumWind = tmpGchNumWindNDR;
						tmpGchCTDepthWind = tmpGchCTDepthWindNDR;
						tmpGchDotWind = tmpGchDotWindNDR;
						tmpGchDepthWind = tmpGchDepthWindNDR;
						tmpGchMethyWind = tmpGchMethyWindNDR;
					}
				}
				else{
					//if(objPre.numValidGch >= NAC.minGchNumLinkerWindow){
					if(objPre.numGchDot >= NAC.minGchDotLinkerWindow && objPre.numValidGch >= NAC.minGchNumLinkerWindow){
						double sigValue = getSigTest(objMid.dot, objPre.dot);
						if(sigValue < NAC.sigValue){
							NDRStartflag = true;
							NDRInflag = true;
							NDRPreLinkerflag = true;
							sigValueMem = sigValue;
							NDRStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
							tmpGchNumWindNDR = objMid.numValidGch;
							tmpGchCTDepthWindNDR = objMid.sumGchCTDepth;
							tmpGchDotWindNDR = objMid.numGchDot;
							tmpGchDepthWindNDR = objMid.sumGchDepth;
							tmpGchMethyWindNDR = objMid.sumGchMethy;
							if(NAC.vm){
								Iterator<NDRCallContext> it = windows.windowsMid.iterator();
								while(it.hasNext()){
								//	if(it.next().getCytosinePatternFlag()){
										Double[] tmpGch = new Double[2];
										NDRCallContext tmp = it.next();
										tmpGch[0] = (double) tmp.getLoc().getStart();
										int[] numVm = validateGch(tmp);
										
										if(numVm != null){
											tmpGch[1] =	(double)numVm[0]/(double)(numVm[0] + numVm[1]);
											gchListForVerboseMode.add(tmpGch);
										}
										
										
								//	}
									
								}
								
							}
							tmpGchNumWindLinker = objPre.numValidGch;
							tmpGchCTDepthWindLinker = objPre.sumGchCTDepth;
							tmpGchDotWindLinker = objPre.numGchDot;
							tmpGchDepthWindLinker = objPre.sumGchDepth;
							tmpGchMethyWindLinker = objPre.sumGchMethy;
							//force wind start at NDR start place
							winStartflag = true;
							winInflag = true;
							winStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
							tmpGchNumWind = objMid.numValidGch;
							tmpGchCTDepthWind = objMid.sumGchCTDepth;
							tmpGchDotWind = objMid.numGchDot;
							tmpGchDepthWind = objMid.sumGchDepth;
							tmpGchMethyWind = objMid.sumGchMethy;

						}
						if(winForceEndflag){
							NDREndCor = windows.windowsMid.peekLast().getLoc().getStart();
							NDREndflag = true;
						}
						
					}
					//if(objPost.numValidGch >= NAC.minGchNumLinkerWindow){
					if(objPost.numGchDot >= NAC.minGchDotLinkerWindow && objPost.numValidGch >= NAC.minGchNumLinkerWindow){
						double sigValue = getSigTest(objMid.dot, objPost.dot);
						if(sigValue < NAC.sigValue){
							NDREndflag = true;
							NDRPostLinkerflag = true;
							NDREndCor = windows.windowsMid.peekLast().getLoc().getStart();
							if(NDRPreLinkerflag){
								tmpGchNumWindLinker += objPost.numValidGch;
								tmpGchCTDepthWindLinker += objPost.sumGchCTDepth;
								tmpGchDotWindLinker += objPost.numGchDot;
								tmpGchDepthWindLinker += objPost.sumGchDepth;
								tmpGchMethyWindLinker += objPost.sumGchMethy;
							}
							else{
								NDRStartflag = true;
								NDRInflag = true;
								sigValueMem = sigValue;
								NDRStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
								tmpGchNumWindNDR = objMid.numValidGch;
								tmpGchCTDepthWindNDR = objMid.sumGchCTDepth;
								tmpGchDotWindNDR = objMid.numGchDot;
								tmpGchDepthWindNDR = objMid.sumGchDepth;
								tmpGchMethyWindNDR = objMid.sumGchMethy;
								
								if(NAC.vm){
									Iterator<NDRCallContext> it = windows.windowsMid.iterator();
									while(it.hasNext()){
										//if(it.next().getCytosinePatternFlag()){
										Double[] tmpGch = new Double[2];
										NDRCallContext tmp = it.next();
										tmpGch[0] = (double) tmp.getLoc().getStart();
										int[] numVm = validateGch(tmp);
										
										if(numVm != null){
											tmpGch[1] =	(double)numVm[0]/(double)(numVm[0] + numVm[1]);
											gchListForVerboseMode.add(tmpGch);
										}
											
											
										//}
										
									}
									
								}
								tmpGchNumWindLinker = objPost.numValidGch;
								tmpGchCTDepthWindLinker = objPost.sumGchCTDepth;
								tmpGchDotWindLinker = objPost.numGchDot;
								tmpGchDepthWindLinker = objPost.sumGchDepth;
								tmpGchMethyWindLinker = objPost.sumGchMethy;
								
								winStartflag = true;
								winInflag = true;
								winStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
								tmpGchNumWind = objMid.numValidGch;
								tmpGchCTDepthWind = objMid.sumGchCTDepth;
								tmpGchDotWind = objMid.numGchDot;
								tmpGchDepthWind = objMid.sumGchDepth;
								tmpGchMethyWind = objMid.sumGchMethy;
							}
							winEndflag = true;
							winEndCor = NDREndCor;
						}
						if(winForceEndflag){
							winEndCor = windows.windowsMid.peekLast().getLoc().getStart();
							winEndflag = true;
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
			if(NDRStartflag && NDRInflag){
				NDREndflag = true;
			//	NDRStartflag = false;
			//	NDRInflag = false;
				NDREndCor = windows.windowsMid.peekLast().getLoc().getStart();
			}
			
		}
		
		
	}

	public void getValueFromWindow(LinkedList<NDRCallContext> windows){
		double averageGchMethy = Double.NaN;
		double sumGchMethy = 0;
		int numValidGch = 0;
		int sumGchCTDepth = 0;
		Iterator<NDRCallContext> itContext = windows.iterator();
		while(itContext.hasNext()){
			NDRCallContext tmpContext = itContext.next();
			if(!tmpContext.getCytosinePatternFlag())
				continue;
			int numC = 0;
			int numT = 0;
			for( PileupElement p : tmpContext.getRealContext().getBasePileup()){
				boolean negStrand = p.getRead().getReadNegativeStrandFlag();
				if(p.getRead().getDuplicateReadFlag()){ //get rid of duplicate reads
                	continue;
                }
				if(((GATKSAMRecord)p.getRead()).isGoodBase(p.getOffset())){
					if(negStrand){
						if(p.getBase()==BaseUtils.G){
							numC++;
						}
						else if(p.getBase()==BaseUtils.A){
							numT++;
						}
						else{
							
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
							
						}
					}
				}
				//if ( BisSNPUtils.usableBase(p, true) ){
					//sumGchSeqDepth++;
				//	if(p.getBase() == BaseUtils.C)
				//		numC++;
				//	if(p.getBase() == BaseUtils.T)
				//		numT++;
		        //}
			}
			if((numC + numT) >= NAC.minCTDepth){
				numValidGch++;
				sumGchCTDepth += (numC + numT);
				sumGchMethy += (double)numC/(double)(numC + numT);
			//	System.out.println("loc: " + tmpContext.getLoc().getStart() + "\tGchMethy: " + (double)numC/(double)(numC + numT) + "\tnumC: " + numC + "\tnumT: " + numT);
			}
		}
		
		if(numValidGch >= NAC.minGchNum){
			summary.sumGchNumInWindCalledConfidently += numValidGch;
			summary.sumGchCTDepthWind += sumGchCTDepth/(double)numValidGch;
			averageGchMethy = sumGchMethy/(double)numValidGch;
		}
		
		if(!Double.isNaN(averageGchMethy)){
			summary.nWindowsCallable++;
			summary.sumGchMethyWindowsCalledConfidently += averageGchMethy;
			windows.getLast().setGchMethyInWindow(averageGchMethy);
			if(averageGchMethy >= NAC.ndrThreshold){
				summary.nNDRWindowsCalledConfidently++;
				summary.sumGchMethyNDRWindowsCalledConfidently += averageGchMethy;
				summary.sumGchNumInNDRWind += numValidGch;
				summary.sumGchCTDepthInNDRWind += sumGchCTDepth/(double)numValidGch;
			//	writer.add(windows.getLast().getLoc(), averageGchMethy);
			}
			
			System.out.println(windows.getFirst().getLoc().getStart() + "\t" + averageGchMethy);
		}
			
		//System.out.println("numValidGch: " + numValidGch + "\tsumGchMethy: " + sumGchMethy + "\taverageGchMethy: " + averageGchMethy);
		
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
		int sumGchCTDepth = 0;
		int sumGchDepth = 0;
		//ArrayList<Double> dot = new ArrayList<Double>();
		ArrayList<Double> dot = new ArrayList<Double>();
		Iterator<NDRCallContext> itContext = window.iterator();
		while(itContext.hasNext()){
			NDRCallContext tmpContext = itContext.next();
			int[] num = validateGch(tmpContext);
			
			if(num != null){
				int numC = num[0];
				int numT = num[1];
				int numO = num[2];
				numValidGch++;
				sumGchCTDepth += (numC + numT);
				numGchDot += (numC + numT);
				sumGchDepth += (numC + numT + numO);
				sumGchMethy += (double)numC/(double)(numC + numT);
				//for(int i=numC; i>=0;i--){
				//	dot.add(1);
				//}
				//for(int i=numT; i>=0;i--){
				//	dot.add(0);
				//}
				dot.add((double)numC/(double)(numC + numT));
			 }	
			//	System.out.println("loc: " + tmpContext.getLoc().getStart() + "\tGchMethy: " + (double)numC/(double)(numC + numT) + "\tnumC: " + numC + "\tnumT: " + numT);
			
		}
		windowsReturnObject obj = new windowsReturnObject(numValidGch, sumGchDepth, sumGchCTDepth, sumGchMethy, numGchDot, dot);
		return obj;
	}

	public class windowsReturnObject {
		public int numValidGch;
		public int numGchDot;
		public int sumGchCTDepth;
		public int sumGchDepth;
		public double sumGchMethy;
		public ArrayList<Double> dot;
		public windowsReturnObject(int numValidGch,int sumGchDepth, int sumGchCTDepth,double sumGchMethy, int numGchDot, ArrayList<Double> dot){
			this.numValidGch = numValidGch;
			this.sumGchCTDepth = sumGchCTDepth;
			this.sumGchDepth = sumGchDepth;
			this.sumGchMethy = sumGchMethy;
			this.numGchDot = numGchDot;
			this.dot = dot;
		}
	}
	
	public int[] validateGch(NDRCallContext ncc){
		int[] num = new int[3]; //0: number of C; 1: number of T; 2: number of other bases;
		num[0] = 0;
		num[1] = 0;
		num[2] = 0;
		if(!ncc.getCytosinePatternFlag())
			return null;
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
			num[0] = numC;
			num[1] = numT;
			num[2] = numO;
			return num;
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
		if(NAC.test.equalsIgnoreCase("ksTest")){
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
	
	
	
	public void addNDRtoWriter(){
		boolean outputFlag = false;
		//if(NDREndflag && tmpGchNumWindNDR  >= NAC.minGchNum && tmpGchCTDepthWindNDR/tmpGchNumWindNDR >= NAC.minCTDepth && tmpGchNumWindLinker  >= NAC.minGchNumLinkerWindow && tmpGchCTDepthWindLinker/tmpGchNumWindLinker >= NAC.minCTDepthLinkerWindow){
		if(NDREndflag && NDRInflag && NDRStartflag){
			ArrayList<Object> valuesNDR = new ArrayList<Object>();
			
			valuesNDR.add(tmpGchMethyWindNDR/tmpGchNumWindNDR);
			valuesNDR.add(tmpGchNumWindNDR);
			valuesNDR.add(tmpGchDepthWindNDR/tmpGchNumWindNDR);
			valuesNDR.add(tmpGchCTDepthWindNDR/tmpGchNumWindNDR);
			valuesNDR.add(tmpGchDotWindNDR);
			valuesNDR.add(tmpGchMethyWindLinker/tmpGchNumWindLinker);
			valuesNDR.add(tmpGchNumWindLinker);
			valuesNDR.add(tmpGchDepthWindLinker/tmpGchNumWindLinker);
			valuesNDR.add(tmpGchCTDepthWindLinker/tmpGchNumWindLinker);
			valuesNDR.add(tmpGchDotWindLinker);
			valuesNDR.add(this.sigValueMem);
			if(NAC.vm){
				Iterator<Double[]> it = gchListForVerboseMode.iterator();
				while(it.hasNext()){
					Double[] tmp = it.next();
				//	System.err.println(tmp[0] + "\t" + tmp[1]);
					valuesNDR.add(tmp[0]);
					valuesNDR.add(tmp[1]);
				}
				
			}
			bedObject bedLineNDR = new bedObject(winChr, NDRStartCor, NDREndCor, valuesNDR); //chr, start, end, aveMethyNDR, gchNumNDR, gchDepthNDR, gchCTdepthNDR, aveMethyLinker, gchNumLinker, gchDepthLinker, gchCTdepthLinker
			writer.add(bedLineNDR);
			
			
			summary.nNDRWindowsCalledConfidently++;
			summary.nNDRWindowsCalledConfidentlyLinker++;
			summary.sumGchCTDepthInNDRWind += tmpGchCTDepthWindNDR;
			summary.sumGchDepthInNDRWind += tmpGchDepthWindNDR;
			summary.sumGchDotInNDRWind += tmpGchDotWindNDR;
			summary.sumGchCTDepthInNDRWindLinker += tmpGchCTDepthWindLinker;
			summary.sumGchDepthInNDRWindLinker += tmpGchDepthWindLinker;
			summary.sumGchDotInNDRWindLinker += tmpGchDotWindLinker;
			summary.sumGchMethyNDRWindowsCalledConfidently += tmpGchMethyWindNDR;
			summary.sumGchMethyNDRWindowsCalledConfidentlyLinker += tmpGchMethyWindLinker;
			summary.sumGchNumInNDRWind += tmpGchNumWindNDR;
			summary.sumGchNumInNDRWindLinker += tmpGchNumWindLinker;
			outputFlag = true;
			if(NAC.ptMode){
				winForceEndflag = true;
			}
			
		}

		if(NAC.ptMode){
			//if(winEndflag && tmpGchNumWind  >= (NAC.minGchNum) && tmpGchCTDepthWind/tmpGchNumWind >= NAC.minCTDepth && this.sigValueMem != -1){
			if(winEndflag && winStartflag && winInflag && winForceEndflag){
				ArrayList<Object> valuesWind = new ArrayList<Object>();
				
				valuesWind.add(tmpGchMethyWind/tmpGchNumWind);
				valuesWind.add(tmpGchNumWind);
				valuesWind.add(tmpGchDepthWind/tmpGchNumWind);
				valuesWind.add(tmpGchCTDepthWind/tmpGchNumWind);
				valuesWind.add(tmpGchDotWind);
				valuesWind.add(this.sigValueMem);
				
				bedObject bedLineWind = new bedObject(winChr, winStartCor, winEndCor, valuesWind); //chr, start, end, aveMethyWind, gchNumWind, gchCTdepthWind
				callableWindWriter.add(bedLineWind);
				summary.nWindowsCallable++;
				summary.sumGchDepthWind += tmpGchDepthWind;
				summary.sumGchCTDepthWind += tmpGchCTDepthWind;
				summary.sumGchMethyWindowsCalledConfidently += tmpGchMethyWind;
				summary.sumGchNumInWindCalledConfidently += tmpGchNumWind;
				summary.sumGchDotWind += tmpGchDotWind;
				
				clearPtModeStatus();
			}
			
		}
		if(outputFlag){
			clearStatus();
			if(NAC.vm){
				gchListForVerboseMode.clear();
			}
		}

	}
	
	private void clearStatus(){
		windows.windowsPre.clear();
		windows.windowsMid.clear();
		windows.windowsPost.clear();
		clearFlagStatus();
		if(NAC.ptMode){
			clearPtModeStatus();
		}
	}
	
	
	private void clearFlagStatus(){
		NDRInflag = false;
		NDRStartflag = false;
		NDREndflag = false;
		NDRPostLinkerflag = false;
		NDRPreLinkerflag = false;
		NDRStartCor = -1;
		NDREndCor = -1;
		winChr = null;
		tmpGchCTDepthWindNDR = 0;
		tmpGchDepthWindNDR = 0;
		tmpGchNumWindNDR = 0;
		tmpGchMethyWindNDR = 0;
		tmpGchDotWindNDR = 0;
		
		tmpGchCTDepthWindLinker = 0;
		tmpGchDepthWindLinker = 0;
		tmpGchNumWindLinker = 0;
		tmpGchMethyWindLinker = 0;
		tmpGchDotWindLinker = 0;
		
		sigValueMem = -1;
		
	}
	
	private void clearPtModeStatus(){
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
	}
	
	private void modifyIntervals(){
		GATKArgumentCollection arg = getToolkit().getArguments();
		List<String> intervalsFileList = arg.intervals;
		 List<String> newIntervals = new ArrayList<String>();
		 Iterator<String> itIntervalFiles = intervalsFileList.iterator();
		 while(itIntervalFiles.hasNext()){
			 BufferedReader br = null;
			 String fileName = itIntervalFiles.next();
			try {
				br = new BufferedReader(new FileReader(fileName));
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			String line;
			String fn = fileName + ".tmp";
			PrintWriter outWriter = null;
			try {
				outWriter = new PrintWriter(new File(fn));
			} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
				e.printStackTrace();
			}
			try {
				while( (line = br.readLine()) != null){
					 String[] tmpArray = line.split("\t");
					// System.err.println(tmpArray[0] + "\t" + tmpArray[1] + "\t" + tmpArray[2] + "\t" + tmpArray.length);
					 Integer start = Integer.parseInt(tmpArray[1]) - NAC.nucLinkerWindow;
					 Integer end = Integer.parseInt(tmpArray[2]) + NAC.nucLinkerWindow;
					 if(start < 0){
						 continue;
					 }
					 String genomeLoc = tmpArray[0] + "\t" + start.toString() + "\t" + end.toString();
					// newIntervals.add(genomeLoc);
					 outWriter.println(genomeLoc);
				 }
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			 
			try {
				br.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			outWriter.close();
			newIntervals.add(fn);
			 
		 }
		 
		 arg.intervals = newIntervals;
		 getToolkit().setArguments(arg);
		// System.err.println(getToolkit().getArguments().intervals.get(0));
	}
	
	private boolean passWindowFilter(NDRCallContext value){
		
		return true;
	}
	/*
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
	*/
}
