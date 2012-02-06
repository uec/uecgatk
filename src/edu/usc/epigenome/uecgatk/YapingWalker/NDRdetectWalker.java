package edu.usc.epigenome.uecgatk.YapingWalker;

/*
 * program measure NDR from the lef first position with sig value to the right first position with sig value. so NDR position would have a little shift (may be 15bp) to the left position. 
 */

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;


import jsc.independentsamples.MannWhitneyTest;
import jsc.independentsamples.SmirnovTest;
import jsc.tests.H1;

import net.sf.samtools.SAMRecord;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisSNPUtils;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisulfiteDiploidSNPGenotypePriors;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisulfiteVariantCallContext;
import org.broadinstitute.sting.gatk.uscec.writer.bedObject;
import org.broadinstitute.sting.gatk.uscec.writer.bedObjectWriterImp;
import org.broadinstitute.sting.gatk.walkers.BAQMode;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Downsample;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
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
import org.broadinstitute.sting.gatk.uscec.YapingWalker.NDRdetectWalker.windowsObject;

@BAQMode(QualityMode = BAQ.QualityMode.OVERWRITE_QUALS, ApplicationTime = BAQ.ApplicationTime.ON_INPUT)
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
	
	private String winChr = null;
	
	private int winStartCor = -1;
			
	private int winEndCor = -1;
	
	private int tmpGchCTDepthWind = 0;
	
	private int tmpGchNumWind = 0;
	

	
	private double tmpGchMethyWind = 0;
	
	private int tmpGchCTDepthWindLinker = 0;
	
	private int tmpGchNumWindLinker = 0;
	
	private double tmpGchMethyWindLinker = 0;
	
	private int tmpGchNumWindNDR = 0;
	
	private int tmpGchCTDepthWindNDR = 0;
	
	private double tmpGchMethyWindNDR = 0;
	
	private double sigValueMem = -1;
	
	public NDRdetectWalker() {
		// TODO Auto-generated constructor stub
		
	}

	 public static class ContextCondition {
		 	long nWindowsVisited = 0;

	       
		 	/** The number of windows called confidently (according to user threshold), contained enough confidantly GCH and enough seq depth and have at least one adjacent window owns enough confidantly GCH and enough seq depth*/
	        long nWindowsCallable = 0;
	        
	        /** The average sequence depth inside window */
	        double sumGchCTDepthWind = 0;
	     
	        /** The number of Gch bases called confidently (according to user threshold), either ref or other */
	        long sumGchNumInWindCalledConfidently = 0;
	        
	        /** The sum of GCH methylation value of windows called confidently (according to user threshold), contained enough confidantly GCH and enough seq depth */
	        double sumGchMethyWindowsCalledConfidently = 0;
	        
	        /** The average sequence depth inside NDR window */
	        double sumGchCTDepthInNDRWind = 0;
	     
	        /** The number of Gch bases called confidently (according to user threshold), either ref or other in NDR windows*/
	        long sumGchNumInNDRWind = 0;
	        
	        /** The number of NDR windows called confidently (according to user threshold), contained enough confidantly GCH and enough seq depth */
	        long nNDRWindowsCalledConfidently = 0;
	        
	        /** The sum of GCH methylation value of windows called confidently (according to user threshold), contained enough confidantly GCH and enough seq depth */
	        double sumGchMethyNDRWindowsCalledConfidently = 0;
	        
	        /** The average sequence depth inside NDR window */
	        double sumGchCTDepthInNDRWindLinker = 0;
	     
	        /** The number of Gch bases called confidently (according to user threshold), either ref or other in NDR windows*/
	        long sumGchNumInNDRWindLinker = 0;
	        
	        /** The number of NDR windows called confidently (according to user threshold), contained enough confidantly GCH and enough seq depth */
	        long nNDRWindowsCalledConfidentlyLinker = 0;
	        
	        /** The sum of GCH methylation value of windows called confidently (according to user threshold), contained enough confidantly GCH and enough seq depth */
	        double sumGchMethyNDRWindowsCalledConfidentlyLinker = 0;
	        
	        
	        
	        double percentCallableWindowOfAll() { return (double)nWindowsCallable/nWindowsVisited;}
	        double percentGchMethyOfCallableWindows() { return (double)sumGchMethyWindowsCalledConfidently/(nWindowsCallable * percentGchNumOfCallableWindows());}
	        double percentGchCTDepthOfCallableWindows() { return (double)sumGchCTDepthWind/(nWindowsCallable * percentGchNumOfCallableWindows());}
	        double percentGchNumOfCallableWindows() { return (double)sumGchNumInWindCalledConfidently/(nWindowsCallable);}
	        double percentNDRWindowOfCallableWindows() { return (double)nNDRWindowsCalledConfidently/nWindowsCallable;}
	        double percentGchMethyOfNDRWindowsCalledConfidently() { return (double)sumGchMethyNDRWindowsCalledConfidently/(nNDRWindowsCalledConfidently * percentGchNumOfNDRWindows());}
	        double percentGchCTDepthOfNDRWindows() { return (double)sumGchCTDepthInNDRWind/(nNDRWindowsCalledConfidently * percentGchNumOfNDRWindows());}
	        double percentGchNumOfNDRWindows() { return (double)sumGchNumInNDRWind/(nNDRWindowsCalledConfidently);}
	        double percentGchMethyOfCallableWindowsLinker() { return (double)sumGchMethyNDRWindowsCalledConfidentlyLinker/(nNDRWindowsCalledConfidentlyLinker * percentGchNumOfCallableWindowsLinker());}
	        double percentGchCTDepthOfCallableWindowsLinker() { return (double)sumGchCTDepthInNDRWindLinker/(nNDRWindowsCalledConfidentlyLinker * percentGchNumOfCallableWindowsLinker());}
	        double percentGchNumOfCallableWindowsLinker() { return (double)sumGchNumInNDRWindLinker/(nNDRWindowsCalledConfidentlyLinker);}
	 }
	 
	 public void initialize() {
	//	 NDRD_engine = new NDRdetectionEngine(getToolkit(), NAC, logger);
		 genotypePriors = new BisulfiteDiploidSNPGenotypePriors();
		 File fn = new File(NAC.outFile);
		 writer = new bedObjectWriterImp(fn);
		 String bedHeadLine = "chr\tstart\tend\taveMethyNDR\tgchNumNDR\tgchCTdepthNDR\taveMethyLinker\tgchNumLinker\tgchCTdepthLinker\tsigValue\n";
		 writer.addHeader(bedHeadLine);
		 if(NAC.ptMode){
			 File fncw = new File(NAC.ocwf);
			 callableWindWriter = new bedObjectWriterImp(fncw);
			 String bedHeadLineWind = "chr\tstart\tend\taveMethyWind\tgchNumWind\tgchCTdepthWind\tsigValue\n";
			 callableWindWriter.addHeader(bedHeadLineWind);
		 }
		 summary = new ContextCondition();
	 }
	
	@Override
	public NDRCallContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
		// TODO Auto-generated method stub
		String cytosinePattern = "GCH-2";
		double methyStatus = 0.453;
		BisSNPUtils it = new BisSNPUtils(NAC);
		AlignmentContext stratifiedContexts = it.getFilteredAndStratifiedContexts(NAC, ref, context);
		
		NDRCallContext ncc = new NDRCallContext(context, ref.getLocus()); // this make some biad for seq poor region's nuc window 
		
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
					windows.windowsPre.offerLast(value);
					
					return windows;
				}
			}
			windows.windowsMid.offerLast(value);
			return windows;
		}
		else if(windows.windowsMid.size() == NAC.nucPosWindow - 1){
			if(windows.windowsMid.peekLast() != null){
				if(value.getLoc().discontinuousP(windows.windowsMid.peekLast().getLoc())){ // there is potential bias here for not enough reads in postWindow..
					clearStatus();
					windows.windowsPre.offerLast(value);
					
					return windows;
				}
			}
			windows.windowsMid.offerLast(value);
		}
		else if(windows.windowsPost.size() <= NAC.nucLinkerWindow - 2){
			if(windows.windowsPost.peekLast() != null){
				if(value.getLoc().discontinuousP(windows.windowsPost.peekLast().getLoc())){
					
					
					winForceEndflag = true;
					if(NAC.statTest){
						getNDRFromWindowsBySigTest(windows);
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
					
					winForceEndflag = true;
					if(NAC.statTest){
						getNDRFromWindowsBySigTest(windows);
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
						getNDRFromWindowsBySigTest(windows);
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
			
		}
		//summary.nWindowsVisited++;
		if(NAC.statTest){
			getNDRFromWindowsBySigTest(windows);
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
		logger.info(String.format("Average GCH CT reads depth in callable windows                                %.2f", summary.percentGchCTDepthOfCallableWindows()));
		logger.info(String.format("Average GCH number in callable windows                                %.2f", summary.percentGchNumOfCallableWindows()));
		
		logger.info(String.format("Average GCH methylation in NDR windows                                %.2f", summary.percentGchMethyOfNDRWindowsCalledConfidently()));
		logger.info(String.format("Average GCH CT reads depth in NDR windows                                %.2f", summary.percentGchCTDepthOfNDRWindows()));
		logger.info(String.format("Average GCH number in NDR windows                                %.2f", summary.percentGchNumOfNDRWindows()));
		
		logger.info(String.format("Average GCH methylation in NDR linker windows                                %.2f", summary.percentGchMethyOfCallableWindowsLinker()));
		logger.info(String.format("Average GCH CT reads depth in NDR linker windows                                %.2f", summary.percentGchCTDepthOfCallableWindowsLinker()));
		logger.info(String.format("Average GCH number in NDR linker windows                                %.2f", summary.percentGchNumOfCallableWindowsLinker()));
		
		writer.close();
		if(NAC.ptMode){
			callableWindWriter.close();
		}
	}
	// by statistics test
	public void getNDRFromWindowsBySigTest(windowsObject windows){
		
		
		windowsReturnObject objMid = getGchListFromWindow(windows.windowsMid);	
		windowsReturnObject objPre = getGchListFromWindow(windows.windowsPre);		
		windowsReturnObject objPost = getGchListFromWindow(windows.windowsPost);
		
		
		
		
		//something wrong here, callable window are not so many. 
		
		if(objMid.numValidGch >= NAC.minGchNum){
			
			
			
			int[] num = validateGch(windows.windowsMid.peekLast());
			int numC = 0;
			int numT = 0;
			if(num != null){
				numC = num[0];
				numT = num[1];
			}
			winChr = windows.windowsMid.peekFirst().getLoc().getContig();
			//look at windows callable or not
			if(!winEndflag){
				if(winStartflag && winInflag){
					if(num != null){
						tmpGchNumWind++;
						tmpGchCTDepthWind = tmpGchCTDepthWind + numC + numT;
						tmpGchMethyWind += (double)numC/(double)(numC + numT);
					}
					if(winForceEndflag){
						winEndflag = true;
						winEndCor = windows.windowsMid.peekLast().getLoc().getStart();
					}
				}
				else{
					if(objPre.numValidGch >= NAC.minGchNumLinkerWindow){
						winPreLinkerflag = true;
						winStartflag = true;
						winInflag = true;
						winStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
						tmpGchNumWind = objMid.numValidGch;
						tmpGchCTDepthWind = objMid.sumGchCTDepth;
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
					if(objPost.numValidGch >= NAC.minGchNumLinkerWindow){
						winPostLinkerflag = true;
						
						if(!winPreLinkerflag){
							winStartflag = true;
							winInflag = true;
							winStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
							tmpGchNumWind = objMid.numValidGch;
							tmpGchCTDepthWind = objMid.sumGchCTDepth;
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
			if(NDRStartflag && NDRInflag){
				if(num != null){
					tmpGchNumWindNDR++;
					tmpGchCTDepthWindNDR = tmpGchCTDepthWindNDR + numC + numT;
					tmpGchMethyWindNDR += (double)numC/(double)(numC + numT);
				}
				if(objPost.numValidGch >= NAC.minGchNumLinkerWindow){
					double sigValue = getSigTest(objMid.dot, objPost.dot);
					if(sigValue < NAC.sigValue){
						NDREndflag = true;
						NDRPostLinkerflag = true;
						NDREndCor = windows.windowsMid.peekLast().getLoc().getStart();
						if(NDRPreLinkerflag){
							tmpGchNumWindLinker += objPost.numValidGch;
							tmpGchCTDepthWindLinker += objPost.sumGchCTDepth;
							tmpGchMethyWindLinker += objPost.sumGchMethy;
						}
						else{
							
							sigValueMem = sigValue;

							tmpGchNumWindLinker = objPost.numValidGch;
							tmpGchCTDepthWindLinker = objPost.sumGchCTDepth;
							tmpGchMethyWindLinker = objPost.sumGchMethy;
							
						}	
						winStartflag = true;
						winInflag = true;
						winEndflag = true;
						winEndCor = NDREndCor;
						winStartCor = NDRStartCor;
						tmpGchNumWind = tmpGchNumWindNDR;
						tmpGchCTDepthWind = tmpGchCTDepthWindNDR;
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
					tmpGchMethyWind = tmpGchMethyWindNDR;
				}
			}
			else{
				if(objPre.numValidGch >= NAC.minGchNumLinkerWindow){
					double sigValue = getSigTest(objMid.dot, objPre.dot);
					if(sigValue < NAC.sigValue){
						NDRStartflag = true;
						NDRInflag = true;
						NDRPreLinkerflag = true;
						sigValueMem = sigValue;
						NDRStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
						tmpGchNumWindNDR = objMid.numValidGch;
						tmpGchCTDepthWindNDR = objMid.sumGchCTDepth;
						tmpGchMethyWindNDR = objMid.sumGchMethy;
						tmpGchNumWindLinker = objPre.numValidGch;
						tmpGchCTDepthWindLinker = objPre.sumGchCTDepth;
						tmpGchMethyWindLinker = objPre.sumGchMethy;
						//force wind start at NDR start place
						winStartflag = true;
						winInflag = true;
						winStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
						tmpGchNumWind = objMid.numValidGch;
						tmpGchCTDepthWind = objMid.sumGchCTDepth;
						tmpGchMethyWind = objMid.sumGchMethy;

					}
					if(winForceEndflag){
						NDREndCor = windows.windowsMid.peekLast().getLoc().getStart();
						NDREndflag = true;
					}
					
				}
				if(objPost.numValidGch >= NAC.minGchNumLinkerWindow){
					double sigValue = getSigTest(objMid.dot, objPost.dot);
					if(sigValue < NAC.sigValue){
						NDREndflag = true;
						NDRPostLinkerflag = true;
						NDREndCor = windows.windowsMid.peekLast().getLoc().getStart();
						if(NDRPreLinkerflag){
							tmpGchNumWindLinker += objPost.numValidGch;
							tmpGchCTDepthWindLinker += objPost.sumGchCTDepth;
							tmpGchMethyWindLinker += objPost.sumGchMethy;
						}
						else{
							NDRStartflag = true;
							NDRInflag = true;
							sigValueMem = sigValue;
							NDRStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
							tmpGchNumWindNDR = objMid.numValidGch;
							tmpGchCTDepthWindNDR = objMid.sumGchCTDepth;
							tmpGchMethyWindNDR = objMid.sumGchMethy;
							tmpGchNumWindLinker = objPost.numValidGch;
							tmpGchCTDepthWindLinker = objPost.sumGchCTDepth;
							tmpGchMethyWindLinker = objPost.sumGchMethy;
							
							winStartflag = true;
							winInflag = true;
							winStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
							tmpGchNumWind = objMid.numValidGch;
							tmpGchCTDepthWind = objMid.sumGchCTDepth;
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
			
			/*
			if(NAC.ptMode){
				winStartflag = true;
				if(!winEndflag){
					if(!winInflag){
						winInflag = true;
						winStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
						tmpGchNumWind = objMid.numValidGch;
						tmpGchCTDepthWind = objMid.sumGchCTDepth;
						tmpGchMethyWind = objMid.sumGchMethy;
					}
					else{
						if(num != null){
							tmpGchNumWind++;
							tmpGchCTDepthWind = tmpGchCTDepthWind + numC + numT;
							tmpGchMethyWind += (double)numC/(double)(numC + numT);
						}
						
					}
				}
			}
			*/
			
			/*
			if(!winStartflag && winStartCor == -1){  //temporaily record start position, this maybe the win start position
				winStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
			}
			
			if(!NDRStartflag && NDRStartCor == -1){  //temporaily record start position, this maybe the NDR start position
				NDRStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
			}
			
			if(objPre.numValidGch >= NAC.minGchNumLinkerWindow){
				double sigValue = getSigTest(objMid.dot, objPre.dot);
			//	if(this.sigValueMem == -1){
					if(!winInflag){
						winInflag = true;
						winStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
						winStartflag = true;
						tmpGchNumWind = objMid.numValidGch;
						tmpGchCTDepthWind = objMid.sumGchCTDepth;
						tmpGchMethyWind = objMid.sumGchMethy;
					}
					else{
						if(num != null){
							tmpGchNumWind++;
							tmpGchCTDepthWind = tmpGchCTDepthWind + numC + numT;
							tmpGchMethyWind += (double)numC/(double)(numC + numT);
						}
					}
					winPreLinkerflag = true;
					if(!NDRInflag && winInflag){
						if(this.sigValueMem >= sigValue || this.sigValueMem == -1){
							this.sigValueMem = sigValue;
						}
							
					}
					
			//	}
					
				if(sigValue < NAC.sigValue){
					NDRStartflag = true;
					NDRPreLinkerflag = true;
					sigValueMem = sigValue;
					if(!NDRInflag){  //start track NDR in the first significant position
						NDRInflag = true;
						winInflag = true;
						NDRStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
						winStartCor = NDRStartCor;
						
						tmpGchNumWindLinker = objPre.numValidGch;
						tmpGchCTDepthWindLinker = objPre.sumGchCTDepth;
						tmpGchMethyWindLinker = objPre.sumGchMethy;
							
						tmpGchNumWindNDR = objMid.numValidGch;
						tmpGchCTDepthWindNDR = objMid.sumGchCTDepth;
						tmpGchMethyWindNDR = objMid.sumGchMethy;
							
							
						tmpGchNumWind = objMid.numValidGch;
						tmpGchCTDepthWind = objMid.sumGchCTDepth;
						tmpGchMethyWind = objMid.sumGchMethy;
								
					}
					else{
						if(num != null){
							tmpGchNumWindNDR++;
							tmpGchCTDepthWindNDR = tmpGchCTDepthWindNDR + numC + numT;
							tmpGchMethyWindNDR += (double)numC/(double)(numC + numT);
						}
					}
					
			//		System.out.println("pre: " + windows.windowsMid.getFirst().getLoc().getStart() + "\t" + averageGchMethy + "\t" + sigValue + "\t" + objPre.sumGchMethy/(double)objPre.numValidGch + "\t" + objPre.numValidGch+ "\t" + objPre.sumGchCTDepth/(double)objPre.numValidGch);
				}
				
				
			}
			
			
			
			if(objPost.numValidGch >= NAC.minGchNumLinkerWindow){
				double sigValue = getSigTest(objMid.dot, objPost.dot);
			//	if(this.sigValueMem == -1){
					if(!winInflag){
						winInflag = true;
						winStartflag = true;
						winStartCor = windows.windowsMid.peekFirst().getLoc().getStart();
						tmpGchNumWind = objMid.numValidGch;
						tmpGchCTDepthWind = objMid.sumGchCTDepth;
						tmpGchMethyWind = objMid.sumGchMethy;
					}
					else{
						if(num != null && !winPreLinkerflag){
							tmpGchNumWind++;
							tmpGchCTDepthWind = tmpGchCTDepthWind + numC + numT;
							tmpGchMethyWind += (double)numC/(double)(numC + numT);
						}
					}
					
					
					if(!NDRInflag && winInflag){
						if(this.sigValueMem >= sigValue || this.sigValueMem == -1){
							this.sigValueMem = sigValue;
						}
							
					}
					
			//	}
				if(sigValue < NAC.sigValue && !NDREndflag){ //end track NDR in the right first significant position
					NDREndflag = true;
					winEndflag = true;
					if(!NDRPreLinkerflag){
						this.sigValueMem = sigValue;
					}
					if(winEndflag && winStartflag){ // window started, but the end, the window end is still not determined
						
					}
					winEndCor = windows.windowsMid.peekLast().getLoc().getStart();
					NDREndCor = windows.windowsMid.peekLast().getLoc().getStart();
					if(NDRPreLinkerflag && num != null){
						
						tmpGchNumWindNDR++;
						tmpGchCTDepthWindNDR = tmpGchCTDepthWindNDR + numC + numT;
						tmpGchMethyWindNDR += (double)numC/(double)(numC + numT);
					}
					tmpGchNumWindLinker += objPost.numValidGch;
					tmpGchCTDepthWindLinker += objPost.sumGchCTDepth;
					tmpGchMethyWindLinker += objPost.sumGchMethy;
					//if(NDRPreLinkerflag){
					//	tmpGchNumWindLinker /= 2;
					///	tmpGchCTDepthWindLinker /= 2;
					//	tmpGchMethyWindLinker /= 2;
					//}
					
				//	writer.add(windows.windowsMid.getLast().getLoc(), averageGchMethy);
			//		System.out.println("post: " + windows.windowsMid.getFirst().getLoc().getStart() + "\t" + averageGchMethy + "\t" + sigValue + "\t" + objPost.sumGchMethy/(double)objPost.numValidGch+ "\t" + objPost.numValidGch+ "\t" + objPost.sumGchCTDepth/(double)objPost.numValidGch);
				}

			}
			
			if(NDREndflag && !NDRStartflag){  //NDR stopped, without knowing which is start position
				tmpGchNumWindNDR = tmpGchNumWind;
				tmpGchCTDepthWindNDR = tmpGchCTDepthWind;
				tmpGchMethyWindNDR = tmpGchMethyWind;
			}

			if(NDRStartflag && !NDREndflag && winEndflag){ //NDR started before, but to the window end, still not know the NDR end position.
				NDREndflag = true;
				NDREndCor = windows.windowsMid.peekLast().getLoc().getStart();
			}
			*/
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
		int sumGchCTDepth = 0;
		ArrayList<Double> dot = new ArrayList<Double>();
		Iterator<NDRCallContext> itContext = window.iterator();
		while(itContext.hasNext()){
			NDRCallContext tmpContext = itContext.next();
			int[] num = validateGch(tmpContext);
			
			if(num != null){
				int numC = num[0];
				int numT = num[1];
				numValidGch++;
				sumGchCTDepth += (numC + numT);
				sumGchMethy += (double)numC/(double)(numC + numT);
				dot.add((double)numC/(double)(numC + numT));
			 }	
			//	System.out.println("loc: " + tmpContext.getLoc().getStart() + "\tGchMethy: " + (double)numC/(double)(numC + numT) + "\tnumC: " + numC + "\tnumT: " + numT);
			
		}
		windowsReturnObject obj = new windowsReturnObject(numValidGch, sumGchCTDepth, sumGchMethy, dot);
		return obj;
	}

	public class windowsReturnObject {
		public int numValidGch;
		public int sumGchCTDepth;
		public double sumGchMethy;
		public ArrayList<Double> dot;
		public windowsReturnObject(int numValidGch,int sumGchCTDepth,double sumGchMethy, ArrayList<Double> dot){
			this.numValidGch = numValidGch;
			this.sumGchCTDepth = sumGchCTDepth;
			this.sumGchMethy = sumGchMethy;
			this.dot = dot;
		}
	}
	
	public int[] validateGch(NDRCallContext ncc){
		int[] num = new int[2]; //0: number of C; 1: number of T;
		num[0] = 0;
		num[1] = 0;
		if(!ncc.getCytosinePatternFlag())
			return null;
		int numC = 0;
		int numT = 0;
		for( PileupElement p : ncc.getRealContext().getBasePileup()){
			SAMRecord samRecord = p.getRead();
        	int offset = p.getOffset();
        	if(offset < 0)//is deletion
        		continue;
        	if(NAC.pairedEndMode){
        		try {
					samRecord = (SAMRecord) p.getRead().clone();
				} catch (CloneNotSupportedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
            	boolean Paired = samRecord.getReadPairedFlag();	
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
	        	if (Paired  && !NAC.USE_BADLY_MATED_READS && !samRecord.getProperPairFlag())
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
			
		}
		if((numC + numT) >= NAC.minCTDepth){
			num[0] = numC;
			num[1] = numT;
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
		if(NAC.ksTest){
			SmirnovTest test = new SmirnovTest(mid, adj, H1.GREATER_THAN);
			return test.getSP();	
		}
		else{
			MannWhitneyTest test = new MannWhitneyTest(mid, adj, H1.GREATER_THAN);
			return test.getSP();
		}
		
		//
		
	}
	
	public void addNDRtoWriter(){
		boolean outputFlag = false;
		//if(NDREndflag && tmpGchNumWindNDR  >= NAC.minGchNum && tmpGchCTDepthWindNDR/tmpGchNumWindNDR >= NAC.minCTDepth && tmpGchNumWindLinker  >= NAC.minGchNumLinkerWindow && tmpGchCTDepthWindLinker/tmpGchNumWindLinker >= NAC.minCTDepthLinkerWindow){
		if(NDREndflag && NDRInflag && NDRStartflag){
			ArrayList<Object> valuesNDR = new ArrayList<Object>();
			
			valuesNDR.add(tmpGchMethyWindNDR/tmpGchNumWindNDR);
			valuesNDR.add(tmpGchNumWindNDR);
			valuesNDR.add(tmpGchCTDepthWindNDR/tmpGchNumWindNDR);
			valuesNDR.add(tmpGchMethyWindLinker/tmpGchNumWindLinker);
			valuesNDR.add(tmpGchNumWindLinker);
			valuesNDR.add(tmpGchCTDepthWindLinker/tmpGchNumWindLinker);
			valuesNDR.add(this.sigValueMem);

			bedObject bedLineNDR = new bedObject(winChr, NDRStartCor, NDREndCor, valuesNDR); //chr, start, end, aveMethyNDR, gchNumNDR, gchCTdepthNDR, aveMethyLinker, gchNumLinker, gchCTdepthLinker
			writer.add(bedLineNDR);
			
			
			summary.nNDRWindowsCalledConfidently++;
			summary.nNDRWindowsCalledConfidentlyLinker++;
			summary.sumGchCTDepthInNDRWind += tmpGchCTDepthWindNDR;
			summary.sumGchCTDepthInNDRWindLinker += tmpGchCTDepthWindLinker;
			summary.sumGchMethyNDRWindowsCalledConfidently += tmpGchMethyWindNDR;
			summary.sumGchMethyNDRWindowsCalledConfidentlyLinker += tmpGchMethyWindLinker;
			summary.sumGchNumInNDRWind += tmpGchNumWindNDR;
			summary.sumGchNumInNDRWindLinker += tmpGchNumWindLinker;
			outputFlag = true;
		}

		if(NAC.ptMode){
			//if(winEndflag && tmpGchNumWind  >= (NAC.minGchNum) && tmpGchCTDepthWind/tmpGchNumWind >= NAC.minCTDepth && this.sigValueMem != -1){
			if(winEndflag && winStartflag && winInflag && winForceEndflag){
				ArrayList<Object> valuesWind = new ArrayList<Object>();
				
				valuesWind.add(tmpGchMethyWind/tmpGchNumWind);
				valuesWind.add(tmpGchNumWind);
				valuesWind.add(tmpGchCTDepthWind/tmpGchNumWind);
				valuesWind.add(this.sigValueMem);
				
				bedObject bedLineWind = new bedObject(winChr, winStartCor, winEndCor, valuesWind); //chr, start, end, aveMethyWind, gchNumWind, gchCTdepthWind
				callableWindWriter.add(bedLineWind);
				summary.nWindowsCallable++;
				summary.sumGchCTDepthWind += tmpGchCTDepthWind;
				summary.sumGchMethyWindowsCalledConfidently += tmpGchMethyWind;
				summary.sumGchNumInWindCalledConfidently += tmpGchNumWind;
				
				clearPtModeStatus();
			}
			
		}
		if(outputFlag){
			clearStatus();
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
		tmpGchNumWindNDR = 0;
		tmpGchMethyWindNDR = 0;
		
		tmpGchCTDepthWindLinker = 0;
		tmpGchNumWindLinker = 0;
		tmpGchMethyWindLinker = 0;
		
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
		tmpGchNumWind = 0;
		tmpGchMethyWind = 0;
	}
}
