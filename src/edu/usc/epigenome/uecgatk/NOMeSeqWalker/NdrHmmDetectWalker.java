/**
 * it output segement with following column: chr, start, end, strand, score, number of GCH, number of CT reads, segment length, p value. after R deal with it, then
 * 2 more column, FDR, GCH methylation value
 */
package edu.usc.epigenome.uecgatk.NOMeSeqWalker;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Iterator;

import jsc.contingencytables.ChiSquaredTest;
import jsc.contingencytables.ContingencyTable2x2;
import jsc.contingencytables.FishersExactTest;
import jsc.tests.H1;


import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BinomialDistributionImpl;
import org.broad.tribble.annotation.Strand;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import edu.usc.epigenome.genomeLibs.FisherExactTest;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObject;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObjectWriterImp;
import edu.usc.epigenome.uecgatk.distribution.OpdfBeta;
import edu.usc.epigenome.uecgatk.distribution.OpdfBetaBinomialFactory;
import edu.usc.epigenome.uecgatk.distribution.OpdfBetaFactory;
import edu.usc.epigenome.uecgatk.hmm.BbForwardBackwardScaledCalculator;
import edu.usc.epigenome.uecgatk.hmm.BbViterbiCalculator;
import edu.usc.epigenome.uecgatk.hmm.NdrForwardBackwardScaledCalculator;
import edu.usc.epigenome.uecgatk.hmm.ObservationMethy;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteVCFConstants;
import edu.usc.epigenome.uecgatk.NOMeSeqWalker.OpdfBetaWriter;
import edu.usc.epigenome.uecgatk.NOMeSeqWalker.OpdfBetaReader;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.ObservationDiscrete;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.OpdfDiscrete;
import be.ac.ulg.montefiore.run.jahmm.OpdfDiscreteFactory;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianFactory;
import be.ac.ulg.montefiore.run.jahmm.ViterbiCalculator;
import be.ac.ulg.montefiore.run.jahmm.apps.sample.SimpleExample.Packet;
import be.ac.ulg.montefiore.run.jahmm.io.FileFormatException;
import be.ac.ulg.montefiore.run.jahmm.io.OpdfGaussianReader;
import be.ac.ulg.montefiore.run.jahmm.io.OpdfGaussianWriter;
import be.ac.ulg.montefiore.run.jahmm.io.OpdfReader;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchScaledLearner;
import be.ac.ulg.montefiore.run.jahmm.learn.KMeansLearner;
import be.ac.ulg.montefiore.run.jahmm.toolbox.KullbackLeiblerDistanceCalculator;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jun 13, 2012 6:50:13 PM
 * 
 */
@Reference(window = @Window(start = -200, stop = 200))
public class NdrHmmDetectWalker extends RodWalker<NdrHmmDetectWalker.Datapoint, NdrHmmDetectWalker.Datum> implements TreeReducible<NdrHmmDetectWalker.Datum> {

	@Input(fullName = "vcf", shortName = "vcf", doc = "input vcf file for detection", required = true)
	public RodBinding<VariantContext> vcf;
	
	@Output(fullName = "states_sequence", shortName = "result", doc = "write hidden states of a sequence to a file (no training mode)", required = false)
	public String resultFile = null;
	
	@Output(fullName = "segment_file", shortName = "segment", doc = "write segment of NDR region into a file (no training mode)", required = false)
	public String segmentFile = null;
	
	@Output(fullName = "segment_file_2", shortName = "segment2", doc = "write segment of NDR region into a file (no training mode)", required = false)
	public String segmentFile2 = null;
	
	@Output(fullName = "MPR_segment_file", shortName = "mprSeg", doc = "write segment of MPR/NPR region into a file (no training mode)", required = false)
	public String mprSegFile = null;
	
	@Argument(fullName = "min_ct_coverage", shortName = "minCT", doc = "minimum number of CT reads for count methylation level, default: 1", required = false)
	public int minCT = 1;
	
	@Argument(fullName = "only_training_mode", shortName = "onlyTrain", doc = "only enable the training mode and output HMM parameter, default: not enabled", required = false)
	public boolean onlyTrain = false;
	
	@Argument(fullName = "only_decoding_mode", shortName = "onlyDecode", doc = "only enable the decoding step and output segments, default: not enabled", required = false)
	public boolean onlyDecode = false;
	
	@Argument(fullName = "beta_decoding_mode", shortName = "beta", doc = "using beta model in the decoding step and output segments, default: not enabled", required = false)
	public boolean beta = false;
	
	@Argument(fullName = "kmeans_initiate", shortName = "kmeans", doc = "use k means method to initiate the best initial parameters, default: not enabled", required = false)
	public boolean kmeans = false;
	
	@Argument(fullName = "hmm_file", shortName = "hmm", doc = "read/write HMM model from/to a file, default: read HMM parameters from a file", required = true)
	public String hmmFile = null;
	
	@Argument(fullName = "max_gap_size", shortName = "gap", doc = "max gap size, default: 10000", required = false)
	public int gap = 10000;
	
	@Argument(fullName = "min_data_point", shortName = "dataP", doc = "minimum data point, default: 5", required = false)
	public int dataP = 5;
	
	@Argument(fullName = "tolerence_level", shortName = "tol", doc = "tolerence level for the converge, default: 1e-5", required = false)
	public double tol = 1e-5;
	
	@Argument(fullName = "max_iteration", shortName = "iteration", doc = "maximum number of iteration for the converge, default: 20", required = false)
	public int iteration = 20;
	
	@Argument(fullName = "halfLocal_mode", shortName = "halfLocal", doc = "enable the halfLocal mode for decoding step p value calculation, default: not enabled", required = false)
	public boolean halfLocal = false;
	
	@Argument(fullName = "window_size", shortName = "window", doc = "in the halfLocal mode, define the window size to detect basic level (100 means +/-100 bp around center), default: 1000000", required = false)
	public long window = 1000000;
	
	@Argument(fullName = "minCT_in_window", shortName = "ctInWindow", doc = "in the halfLocal mode, minimum of CT reads in the window, default: 100", required = false)
	public int ctInWindow = 100;
	
	@Argument(fullName = "minNumGch_in_window", shortName = "gchInWindow", doc = "in the halfLocal mode, minimum number of GCH in the window, default: 10", required = false)
	public int gchInWindow = 10;
	
	@Argument(fullName = "num_of_states", shortName = "states", doc = "number of states, default: 2", required = false)
	public static int states = 2;
	
	//private static int STATES=states;
	
	//private double TOLERENCE=tol;
	
	//private int MAXIMUM_GAP_SIZE = gap; // maximum gap allowed for split different sequence for training & decoding
	
	private double BOUNDARY_PROBABILITY = 0.90;
	
	//private int MINIMUM_DATA_POINTS=dataP;
	
	private Hmm<ObservationReal> hmm = null;
	
	private bedObjectWriterImp bedWriter = null;
	
	private bedObjectWriterImp segWriter = null;
	
	private bedObjectWriterImp nprSegWriter = null;
	
	private bedObjectWriterImp segWriter2 = null;
	

	
//	public SlidingWindow slidingWindowLeft;
//	public SlidingWindow slidingWindowRight;
	
	public int bpLeftBound = 0;
	public int bpRightBound = 0 ;
	
	public void initialize(){
		if(!onlyTrain){
			bedWriter = new bedObjectWriterImp(new File(resultFile));
			segWriter = new bedObjectWriterImp(new File(segmentFile));
		//	segWriter2 = new bedObjectWriterImp(new File(segmentFile2));
			nprSegWriter = new bedObjectWriterImp(new File(mprSegFile));
			if(halfLocal){
			//	slidingWindowLeft = new SlidingWindow(ctInWindow, gchInWindow, window);
			//	slidingWindowRight = new SlidingWindow(ctInWindow, gchInWindow, window);
			}
		}
			
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.LocusWalker#map(org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker, org.broadinstitute.sting.gatk.contexts.ReferenceContext, org.broadinstitute.sting.gatk.contexts.AlignmentContext)
	 */
	@Override
	public Datapoint map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
		if (tracker == null)
			return null;
		List<VariantContext> vcf_bindings = tracker.getValues(vcf);
		if (!vcf_bindings.isEmpty()) {

			VariantContext vc = vcf_bindings.get(0);
			if(!vc.isFiltered() && vc.hasGenotypes()){
				if(vc.getGenotype(0).hasAttribute(BisulfiteVCFConstants.BEST_C_PATTERN)){
					if(vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".").equalsIgnoreCase("GCH")){
						int numC = vc.getGenotype(0).getAttributeAsInt(BisulfiteVCFConstants.NUMBER_OF_C_KEY, -1);
						int numT = vc.getGenotype(0).getAttributeAsInt(BisulfiteVCFConstants.NUMBER_OF_T_KEY, -1);
						if (numC != -1 && numT != -1 && (numC + numT >= minCT)){
							double methyValue = (double) numC / (double) (numC + numT);
							
							if(numC == 0)
								methyValue += (Math.random())/1000000000;
							else if(numT == 0)
								methyValue -= (Math.random())/1000000000;
							//double methyValue =(double) numC / (double) (numC + numT);
							//System.err.println(ref.getLocus());
							//System.err.println(new ObservationReal(methyValue));
							if(!onlyTrain && halfLocal){
								Datapoint dat = new Datapoint(ref.getLocus(),new ObservationReal(methyValue), (numC + numT), numC);
								return dat;
							}
							Datapoint dat = new Datapoint(ref.getLocus(),new ObservationReal(methyValue), (numC + numT));
							//dat.value.setCoverage((numC + numT));
							return dat;
						}
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
	public Datum reduceInit() {
		Datum tmp = new Datum();
		LinkedList<GenomeLoc> position = new LinkedList<GenomeLoc>();
		LinkedList<ObservationReal> value = new LinkedList<ObservationReal>();
		LinkedList<Integer> numCT = new LinkedList<Integer>();
		LinkedList<Integer> numC = new LinkedList<Integer>();
		
		LinkedList<ObservationMethy> methy = new LinkedList<ObservationMethy>();
		
		tmp.position = new ArrayList<LinkedList<GenomeLoc>>();
		tmp.value = new ArrayList<LinkedList<ObservationReal>>();
		tmp.methy = new ArrayList<LinkedList<ObservationMethy>>();
		tmp.numCT = new ArrayList<LinkedList<Integer>>();
		tmp.numC = new ArrayList<LinkedList<Integer>>();
		//tmp.numCLeftBound = new HashMap<GenomeLoc, Integer>();
	//	tmp.numCtLeftBound = new HashMap<GenomeLoc, Integer>();
	//	tmp.numCRightBound = new HashMap<GenomeLoc, Integer>();
	//	tmp.numCtRightBound = new HashMap<GenomeLoc, Integer>();
		
		tmp.position.add(position);
		tmp.value.add(value);
		tmp.methy.add(methy);
		tmp.numCT.add(numCT);
		tmp.numC.add(numC);
		return tmp;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Datum reduce(Datapoint value, Datum sum) {
		
		if(value == null)
			return sum;
		
		LinkedList<GenomeLoc> position = sum.position.remove(sum.position.size()-1);
		LinkedList<ObservationReal> data = sum.value.remove(sum.value.size()-1);
		LinkedList<ObservationMethy> methy = sum.methy.remove(sum.methy.size()-1);
		LinkedList<Integer> numCT = sum.numCT.remove(sum.numCT.size()-1);
		LinkedList<Integer> numC = null;
		if(!onlyTrain && halfLocal)
			numC = sum.numC.remove(sum.numC.size()-1);
		if(!position.isEmpty()){
		
			//System.err.println(gap + "\t" + position.peekLast().distance(value.position));
				if( position.peekLast().distance(value.position) >= gap || !position.peekLast().onSameContig(value.position)){
					if(position.size() >= dataP){
						sum.position.add(position);
						sum.value.add(data);
						sum.methy.add(methy);
						sum.numCT.add(numCT);
						if(!onlyTrain && halfLocal){
							sum.numC.add(numC);
						}
					}
					
					LinkedList<GenomeLoc> newPosition = new LinkedList<GenomeLoc>();
					LinkedList<ObservationReal> newValue = new LinkedList<ObservationReal>();
					LinkedList<ObservationMethy> newMethy = new LinkedList<ObservationMethy>();
					LinkedList<Integer> newNumCT = new LinkedList<Integer>();
					LinkedList<Integer> newNumC = new LinkedList<Integer>();
					newPosition.offerLast(value.position);
					newValue.offerLast(value.value);
					newNumCT.offerLast(value.numCT);
					newMethy.offerLast(value.methy);
					
					sum.position.add(newPosition);
					sum.value.add(newValue);
					sum.methy.add(newMethy);
					sum.numCT.add(newNumCT);
					if(!onlyTrain && halfLocal){
						newNumC.offerLast(value.numC);
						sum.numC.add(newNumC);
					}
					return sum;
				}
			
			
		}
		position.offerLast(value.position);
		data.offerLast(value.value);
		methy.offerLast(value.methy);
		numCT.offerLast(value.numCT);
		sum.position.add(position);
		sum.value.add(data);
		sum.methy.add(methy);
		sum.numCT.add(numCT);
		if(!onlyTrain && halfLocal){
			numC.offerLast(value.numC);
			sum.numC.add(numC);
		}
		return sum;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.TreeReducible#treeReduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Datum treeReduce(Datum lhs, Datum rhs) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public void onTraversalDone(Datum result) {
		ArrayList<ArrayList<ObservationReal>> listValues = new ArrayList<ArrayList<ObservationReal>>();
		ArrayList<ArrayList<ObservationMethy>> listMethys = new ArrayList<ArrayList<ObservationMethy>>();
		ArrayList<ArrayList<GenomeLoc>> listPositions = new ArrayList<ArrayList<GenomeLoc>>();
		ArrayList<ArrayList<Integer>> listNumCTs = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Integer>> listNumCs = new ArrayList<ArrayList<Integer>>();
		System.out.println("sequence data size: " + result.value.size());
		for(int z = 0; z < result.value.size(); z++){
			System.out.println("size: " + result.value.get(z).size() + "\tlength: " + result.position.get(z).peekLast().distance(result.position.get(z).peekFirst()));
			if(result.value.get(z).size() >= dataP){
			//	result.value.remove(z);
			//	result.methy.remove(z);
			//	result.position.remove(z);
			//	result.numCT.remove(z);
			//	z--;
				listValues.add( new ArrayList<ObservationReal>(result.value.get(z)));
				listMethys.add( new ArrayList<ObservationMethy>(result.methy.get(z)));
				listPositions.add( new ArrayList<GenomeLoc>(result.position.get(z)));
				listNumCTs.add( new ArrayList<Integer>(result.numCT.get(z)));
				if(!onlyTrain && halfLocal){
					listNumCs.add( new ArrayList<Integer>(result.numC.get(z)));
				}
				
			}	
		}
	//	values.add(result.value);
		if(!onlyDecode){
			System.out.println("training....");
			if(kmeans){
				hmm = buildInitHmmByBeta(listValues);
			}
			else{
				hmm = buildInitHmmRandomlyByBeta(listValues);
			}
			
			BaumWelchScaledLearner bwl = new BaumWelchScaledLearner();
//iteration
			Hmm<ObservationReal> prevHmm = null;
			try {
				prevHmm = hmm.clone();
			} catch (CloneNotSupportedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			
			// This object measures the distance between two HMMs
			KullbackLeiblerDistanceCalculator klc = 
				new KullbackLeiblerDistanceCalculator();
			
			double distance = Double.MAX_VALUE;
			int i = 0;
			// Incrementally improve the solution
			while(Math.abs(distance) >= tol){
				i++;
				
				try {
					prevHmm = hmm.clone();
				} catch (CloneNotSupportedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				System.out.println("HMM pre:\n" + prevHmm);
				
				//hmm = bwl.iterate(hmm, result.value);
				hmm = bwl.iterate(hmm, listValues);
				distance = klc.distance(prevHmm, hmm);
				System.out.println("Distance at iteration " + i + ": " +
						distance);
			}  

			System.out.println("Resulting HMM:\n" + hmm);
			try {
				FileWriter writer = new FileWriter(hmmFile);
				HmmWriter.write(writer, new OpdfBetaWriter(), hmm);
				writer.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		else{
			try {
				
				hmm = HmmReader.read(new FileReader(hmmFile), new OpdfBetaReader());
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (FileFormatException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}	
			
		if(!onlyTrain){
				
		
			int ndrState = 0;
		//	int ndrState2 = 1;
		//	double maxMean = Math.max(((OpdfBeta)hmm.getOpdf(0)).mean(), Math.max(((OpdfBeta)hmm.getOpdf(1)).mean(), ((OpdfBeta)hmm.getOpdf(2)).mean()));
			double maxMean = Math.max(((OpdfBeta)hmm.getOpdf(1)).mean(), ((OpdfBeta)hmm.getOpdf(0)).mean());
			int nprState = 2;
		//	double minMean = Math.min(((OpdfBeta)hmm.getOpdf(0)).mean(), Math.min(((OpdfBeta)hmm.getOpdf(1)).mean(), ((OpdfBeta)hmm.getOpdf(2)).mean()));
			double minMean = Math.min(((OpdfBeta)hmm.getOpdf(1)).mean(), ((OpdfBeta)hmm.getOpdf(0)).mean());
			//suppose NDR region have the highest mean GCH methyation value, NPR region have lowest mean GCH methylation value
			for(int m = 0; m < hmm.nbStates(); m++){
				if(((OpdfBeta)hmm.getOpdf(m)).mean() >= maxMean){
					ndrState = m;
				}
				else if(((OpdfBeta)hmm.getOpdf(m)).mean() <= minMean ){

					nprState = m;
				}
				else{
			//		ndrState2 = m;
				}
			}
		//	System.err.println(nprState + "\t" + ndrState2 + "\t" + ndrState);
		//	double[] randomSeqs = getPermutatedSeqs(result.value, hmm, ndrState);
			
			
			for(int j=0; j < listMethys.size(); j++){
			//for(ArrayList<ObservationMethy> methy : listMethys){
			//for(ArrayList<ObservationReal> methy : listValues){
				//int[] hiddenState = hmm.mostLikelyStateSequence(value);
				int[] hiddenState = null;
				if(beta){
					hiddenState = hmm.mostLikelyStateSequence(listValues.get(j));
				}
				else{
					hiddenState = (new BbViterbiCalculator(listMethys.get(j), hmm)).stateSequence();
				}
				

				//int[] hiddenState = (new ViterbiCalculator(methy, hmm)).stateSequence();
				GenomeLoc[] loci = new GenomeLoc[listPositions.get(j).size()];
				Iterator<GenomeLoc> it = listPositions.get(j).iterator();
				int ii=0;
				while(it.hasNext()){
					loci[ii] = it.next();
					ii++;
				}
				
				ii=0;
				double[] methyState = new double[listMethys.get(j).size()];
				Iterator<ObservationMethy> it2 = listMethys.get(j).iterator();
				//Iterator<ObservationReal> it2 = methy.iterator();
				while(it2.hasNext()){
					methyState[ii] = it2.next().value;
					ii++;
				}
				
				ii=0;
				int[] numCTState = new int[listNumCTs.get(j).size()];
				Iterator<Integer> it3 = listNumCTs.get(j).iterator();
				while(it3.hasNext()){
					numCTState[ii] = it3.next();
					ii++;
				}
				
				ii=0;
				int[] numCState = new int[listNumCs.get(j).size()];
				Iterator<Integer> it4 = listNumCs.get(j).iterator();
				while(it4.hasNext()){
					numCState[ii] = it4.next();
					ii++;
				}
				
				
				if(beta){
					NdrForwardBackwardScaledCalculator nfbsc = new NdrForwardBackwardScaledCalculator(listValues.get(j), hmm);
					getNPRSegmentByHalfLocalComparFBSC(hiddenState, methyState, loci, numCTState, numCState, ndrState, segWriter, true, hmm, nfbsc);
					getNPRSegmentByHalfLocalComparFBSC(hiddenState, methyState, loci, numCTState, numCState, nprState, nprSegWriter, false, hmm, nfbsc);
				}
				else{
					BbForwardBackwardScaledCalculator nfbsc = new BbForwardBackwardScaledCalculator(listMethys.get(j), hmm);
					getNPRSegmentByHalfLocalCompar(hiddenState, methyState, loci, numCTState, numCState, ndrState, segWriter, true, hmm, nfbsc);
					getNPRSegmentByHalfLocalCompar(hiddenState, methyState, loci, numCTState, numCState, nprState, nprSegWriter, false, hmm, nfbsc);
				}
				//getNDRSegment(hiddenState, methyState, loci, numCTState, randomSeqs, ndrState);
				//getNPRSegmentByLocalCompar(hiddenState, methyState, loci, numCTState, ndrState, segWriter, true);
				//getNPRSegmentByLocalCompar(hiddenState, methyState, loci, numCTState, nprState, nprSegWriter, false);
				
				//getNPRSegmentByHalfLocalCompar(hiddenState, methyState, loci, numCTState, ndrState, segWriter, true);
				//getNPRSegmentByHalfLocalCompar(hiddenState, methyState, loci, numCTState, nprState, nprSegWriter, false);
				
			//	getNPRSegmentByHalfLocalCompar(hiddenState, methyState, loci, numCTState, numCState, ndrState2, segWriter2, true);
				
				
				
				for(int i = 0; i < hiddenState.length; i++){
					List<Object> tmp = new ArrayList<Object>();
					tmp.add(hiddenState[i]);
					tmp.add(numCState[i]);
					tmp.add(numCTState[i]);
					//tmp.add(result.numCLeftBound.get(loci));
					//tmp.add(result.numCtLeftBound.get(loci));
					//tmp.add(result.numCRightBound.get(loci));
					//tmp.add(result.numCtRightBound.get(loci));
					bedObject bedLine = new bedObject(loci[i].getContig(), loci[i].getStart()-1, loci[i].getStop(), Strand.NONE, (List)tmp);
					bedWriter.add(bedLine);
				}
				//j++;
			}

			bedWriter.close();
			segWriter.close();
		//	segWriter2.close();
			nprSegWriter.close();
		}
	}
	

	
	private static Hmm<ObservationReal> buildInitHmmByBeta(ArrayList<ArrayList<ObservationReal>> seqs)
	{	

		KMeansLearner<ObservationReal> kl = new KMeansLearner<ObservationReal>(states, new OpdfBetaFactory(),
				seqs);
		System.out.println("KMeansLearner...");
		Hmm<ObservationReal> hmm = kl.learn();
		return hmm;
	}

	private static Hmm<ObservationReal> buildInitHmmRandomlyByBeta(ArrayList<ArrayList<ObservationReal>> seqs)
	{	
		Hmm<ObservationReal> hmm = 
				new Hmm<ObservationReal>(2,new OpdfBetaFactory());
		double tmp3 = Math.random();
			hmm.setPi(0, 0.8);
			hmm.setPi(1, 0.2);
			
			hmm.setOpdf(0, new OpdfBeta(0.03,1.1));
			hmm.setOpdf(1, new OpdfBeta(1.0,0.3));
			double tmp1 = Math.random();
			
			hmm.setAij(0, 1, 0.1);
			hmm.setAij(0, 0, 0.9);
			double tmp2 = Math.random();
			hmm.setAij(1, 0, 0.1);
			hmm.setAij(1, 1, 0.9);
		return hmm;

	}
		
		private void getNPRSegmentByHalfLocalCompar(int[] hiddenState, double[] methyState, GenomeLoc[] loci, int[] numCTState, int[] numCState, int nprState, bedObjectWriterImp writer, boolean reverseP, Hmm<ObservationReal> hmm, BbForwardBackwardScaledCalculator nfbsc){
			//hash the ct reads in the adjacent window information
			
		
			HashMap<GenomeLoc, Integer> numCtLeftBound = new HashMap<GenomeLoc, Integer>();
			HashMap<GenomeLoc, Integer> numCLeftBound = new HashMap<GenomeLoc, Integer>();
			HashMap<GenomeLoc, Integer> numCtRightBound = new HashMap<GenomeLoc, Integer>();
			HashMap<GenomeLoc, Integer> numCRightBound = new HashMap<GenomeLoc, Integer>();
			//System.err.println(numCTState.length + "\t" + numCState.length);
			SlidingWindow slidingWindowLeft = new SlidingWindow(ctInWindow, gchInWindow, window, nprState);
			SlidingWindow slidingWindowRight = new SlidingWindow(ctInWindow, gchInWindow, window, nprState);
			
			for(int z=0; z < hiddenState.length; z++){
				Datapoint data = new Datapoint(loci[z],new ObservationReal(methyState[z]),numCTState[z],numCState[z],hiddenState[z]);
			///	if(loci[z].getStart() == 2007766)
			//		System.err.println(numCRightBound.size() + "\t" + slidingWindowRight.getCReadsNum() + "\t" + slidingWindowRight.getLength() + "\t" + slidingWindowLeft.getLength() + "\t" + z);
				if(slidingWindowLeft.getLength() < window || slidingWindowLeft.getGchNum() < gchInWindow || slidingWindowLeft.getCtReadsNum() < ctInWindow){ // in the beginning of the chromosome
					slidingWindowLeft.addLast(data);
					if(z==hiddenState.length-1){
						int numC = slidingWindowLeft.getCReadsNum();
						int numCt = slidingWindowLeft.getCtReadsNum();
						for(int j = 0; j <= z; j++){
							numCLeftBound.put(loci[j], numC);
							numCtLeftBound.put(loci[j], numCt);
							numCRightBound.put(loci[j], 0);
							numCtRightBound.put(loci[j], 0);
						}
					}
					
				//	System.err.println("left: " + slidingWindowLeft.getLength() + "\t" + data.position);
				}
				else if(slidingWindowRight.windowList.isEmpty()){
					int numC = slidingWindowLeft.getCReadsNum();
					int numCt = slidingWindowLeft.getCtReadsNum();
					for(int j = 0; j < z; j++){
						numCLeftBound.put(loci[j], numC);
						numCtLeftBound.put(loci[j], numCt);
						numCRightBound.put(loci[j], 0);
						numCtRightBound.put(loci[j], 0);
					}
					slidingWindowRight.addLast(data);
					numCLeftBound.put(data.position, numC);
					numCtLeftBound.put(data.position, numCt);
					numCRightBound.put(data.position, 0);
					numCtRightBound.put(data.position, 0);
				}
				else if(slidingWindowRight.getLength() < window || slidingWindowRight.getGchNum() < gchInWindow || slidingWindowRight.getCtReadsNum() < ctInWindow){
					slidingWindowRight.addLast(data);
					numCLeftBound.put(data.position, slidingWindowRight.getCReadsNum());
					numCtLeftBound.put(data.position, slidingWindowRight.getCtReadsNum());
					numCRightBound.put(data.position, 0);
					numCtRightBound.put(data.position, 0);
					if(z==hiddenState.length-1){
						int numC = slidingWindowRight.getCReadsNum();
						int numCt = slidingWindowRight.getCtReadsNum();
						for(int j = z - slidingWindowRight.windowList.size(); j <= z; j++){
							numCRightBound.put(loci[j], numC);
							numCtRightBound.put(loci[j], numCt);
						}
					}
					
				//	System.err.println("right: " + slidingWindowRight.getLength() + "\t" + data.position);
				}
				else{
					int numC = slidingWindowRight.getCReadsNum();
					int numCt = slidingWindowRight.getCtReadsNum();
					LinkedList<Datapoint> tmpValue = slidingWindowRight.addLast(data, true);
					numCLeftBound.put(data.position, slidingWindowRight.getCReadsNum());
					numCtLeftBound.put(data.position, slidingWindowRight.getCtReadsNum());
					numCRightBound.put(data.position, 0);
					numCtRightBound.put(data.position, 0);
					for(Datapoint tmpData: tmpValue){
						
						numCLeftBound.put(tmpData.position, slidingWindowLeft.getCReadsNum());
						numCtLeftBound.put(tmpData.position, slidingWindowLeft.getCtReadsNum());
						slidingWindowLeft.addLast(tmpData, true);
						numCRightBound.put(tmpData.position, numC);
						numCtRightBound.put(tmpData.position, numCt);
						//numC -= tmpData.numC;
						//numCt -= tmpData.numCT;
						

						
					}
					// first add CT reads number when new value add into slidingWindowRight, it will be updated when they are popped out.
					
					
				//	System.err.println("left: " + slidingWindowLeft.getLength() + "\t" + data.position + "\t" + slidingWindowLeft.getCReadsNum() + "\t" + slidingWindowLeft.getCtReadsNum() + "\t" + numC);
				//	System.err.println("right: " + slidingWindowRight.getLength() + "\t" + data.position + "\t" + slidingWindowRight.getCReadsNum() + "\t" + slidingWindowRight.getCtReadsNum() + "\t" + numCt);
				}
			}
			for(GenomeLoc loc : numCtRightBound.keySet()){
				Integer ct = numCtRightBound.get(loc);
				Integer c = numCRightBound.get(loc);
				if(ct < c || c <0 || ct <0)
					System.err.println(c + "\t" + ct + "\t" + loc);
			}
			
			String chr = null;
			int start = -1;
			int end = -1;
			double score = 0;
			int dataPoint = 0; //number of GCH in the segments
			int numC_npr = 0;
			int numT_npr = 0;
			//double sumMethy_Gch_back = 0;
			//int numGch_back = 0;
			//int numC_back = 0;
			//int numT_back = 0;
			
			GenomeLoc startLoc = null;
			GenomeLoc preLoc = null;
			int preState = -1; //record the previous loci's state
			
			double preWeight = -1;
			double postWeight = -1;
			
			for(int i=0; i < hiddenState.length; i++){
				if(preState == -1){
					preState = hiddenState[i]; 
					if(preState == nprState){
						chr = loci[i].getContig();
						//start = loci[i].getStart()-1;
						start = loci[i].getStart()-1;
						startLoc = loci[i];
						score = methyState[i];
						dataPoint=1;
						numC_npr = numCState[i];
						numT_npr = numCTState[i] - numCState[i];

					}	
					
				}
				else{
					if(preState != nprState){
						if(hiddenState[i] == nprState){
								chr = loci[i].getContig();
								preWeight = 1 - nfbsc.getAlpha(hmm, i, preState, hiddenState[i]);
								//double weight = hmm.getAij(preState, hiddenState[i]);
								//start =(int)((loci[i].getStart() + preLoc.getStart())/2) ;
								start =(int)( (loci[i].getStart() - preLoc.getStart()) * preWeight + preLoc.getStart()) ;
								//start =preLoc.getStart() ;
								startLoc = loci[i];
								numC_npr = numCState[i];
								numT_npr = numCTState[i] - numCState[i];
								score = methyState[i];
								dataPoint = 1;
						}

					}
					else{
						if(hiddenState[i] == nprState){
							numC_npr += numCState[i];
							numT_npr += (numCTState[i] - numCState[i]);
							score += methyState[i];
							dataPoint++;
						}
						else{
						//	System.err.println(preLoc);
							//end = (int)((loci[i].getStart() + preLoc.getStart())/2);
							//if(i < hiddenState.length -1){
								postWeight = 1 - nfbsc.getAlpha(hmm, i, preState, hiddenState[i]);
							//}
							//else{
							//	postWeight = 1;
							//}
							//double weight = hmm.getAij(preState, hiddenState[i]);
							end =(int)( (loci[i].getStart() - preLoc.getStart()) * postWeight + preLoc.getStart()) ;
							//end =preLoc.getStart();
							//int numGch_back = backGround.numCLeftBound.size();
							
						//	System.err.println(startLoc);
						//	System.err.println(numCLeftBound.size());
						//	System.err.println(numCRightBound.size());
						//	System.err.println(numCLeftBound.get(startLoc));
						//	System.err.println(numCtLeftBound.get(startLoc));
						//	System.err.println(numCRightBound.get(preLoc));
						//	System.err.println(numCtRightBound.get(preLoc)); 
							int numC_back = numCLeftBound.get(startLoc) + numCRightBound.get(preLoc);
							int numCT_back = numCtLeftBound.get(startLoc) + numCtRightBound.get(preLoc);
							double pValue = getBinomialSigTest(numC_npr, (numC_npr + numT_npr), (double)numC_back/(double)numCT_back, reverseP);
							List<Object> tmp = new LinkedList<Object>();
							
							tmp.add(score);
							tmp.add(dataPoint);
							tmp.add(numC_npr);
							tmp.add(numT_npr);
							tmp.add(numC_back);
							tmp.add(numCT_back);
							//tmp.add(sumMethy_Gch_back);
							//tmp.add(numGch_back);
							tmp.add(end-start);
							tmp.add(1-preWeight);
							tmp.add(1-postWeight);
							tmp.add(pValue);
						//	if((numC_npr + numT_npr >= minCT) && (numC_back + numT_back) >= minCT){
							if((numC_npr + numT_npr) >= minCT && numCT_back >= minCT){
								bedObject bedLine = new bedObject(chr, start, end, Strand.NONE, (List)tmp);
								writer.add(bedLine);
							}
						//	}
						}
					}
					preState = hiddenState[i];
				}
				preLoc = loci[i];
			}
				
					
		}
		
		//boundary not by viterbi decoding, but by forward and backward scaled calculator
		private void getNPRSegmentByHalfLocalComparFBSC(int[] hiddenState, double[] methyState, GenomeLoc[] loci, int[] numCTState, int[] numCState, int nprState, bedObjectWriterImp writer, boolean reverseP, Hmm<ObservationReal> hmm, NdrForwardBackwardScaledCalculator nfbsc){
			//hash the ct reads in the adjacent window information
			
		
			HashMap<GenomeLoc, Integer> numCtLeftBound = new HashMap<GenomeLoc, Integer>();
			HashMap<GenomeLoc, Integer> numCLeftBound = new HashMap<GenomeLoc, Integer>();
			HashMap<GenomeLoc, Integer> numCtRightBound = new HashMap<GenomeLoc, Integer>();
			HashMap<GenomeLoc, Integer> numCRightBound = new HashMap<GenomeLoc, Integer>();
			//System.err.println(numCTState.length + "\t" + numCState.length);
			SlidingWindow slidingWindowLeft = new SlidingWindow(ctInWindow, gchInWindow, window, nprState);
			SlidingWindow slidingWindowRight = new SlidingWindow(ctInWindow, gchInWindow, window, nprState);
			
			for(int z=0; z < hiddenState.length; z++){
				Datapoint data = new Datapoint(loci[z],new ObservationReal(methyState[z]),numCTState[z],numCState[z],hiddenState[z]);
			///	if(loci[z].getStart() == 2007766)
			//		System.err.println(numCRightBound.size() + "\t" + slidingWindowRight.getCReadsNum() + "\t" + slidingWindowRight.getLength() + "\t" + slidingWindowLeft.getLength() + "\t" + z);
				if(slidingWindowLeft.getLength() < window || slidingWindowLeft.getGchNum() < gchInWindow || slidingWindowLeft.getCtReadsNum() < ctInWindow){ // in the beginning of the chromosome
					slidingWindowLeft.addLast(data);
					if(z==hiddenState.length-1){
						int numC = slidingWindowLeft.getCReadsNum();
						int numCt = slidingWindowLeft.getCtReadsNum();
						for(int j = 0; j <= z; j++){
							numCLeftBound.put(loci[j], numC);
							numCtLeftBound.put(loci[j], numCt);
							numCRightBound.put(loci[j], 0);
							numCtRightBound.put(loci[j], 0);
						}
					}
					
				//	System.err.println("left: " + slidingWindowLeft.getLength() + "\t" + data.position);
				}
				else if(slidingWindowRight.windowList.isEmpty()){
					int numC = slidingWindowLeft.getCReadsNum();
					int numCt = slidingWindowLeft.getCtReadsNum();
					for(int j = 0; j < z; j++){
						numCLeftBound.put(loci[j], numC);
						numCtLeftBound.put(loci[j], numCt);
						numCRightBound.put(loci[j], 0);
						numCtRightBound.put(loci[j], 0);
					}
					slidingWindowRight.addLast(data);
					numCLeftBound.put(data.position, numC);
					numCtLeftBound.put(data.position, numCt);
					numCRightBound.put(data.position, 0);
					numCtRightBound.put(data.position, 0);
				}
				else if(slidingWindowRight.getLength() < window || slidingWindowRight.getGchNum() < gchInWindow || slidingWindowRight.getCtReadsNum() < ctInWindow){
					slidingWindowRight.addLast(data);
					numCLeftBound.put(data.position, slidingWindowRight.getCReadsNum());
					numCtLeftBound.put(data.position, slidingWindowRight.getCtReadsNum());
					numCRightBound.put(data.position, 0);
					numCtRightBound.put(data.position, 0);
					if(z==hiddenState.length-1){
						int numC = slidingWindowRight.getCReadsNum();
						int numCt = slidingWindowRight.getCtReadsNum();
						for(int j = z - slidingWindowRight.windowList.size(); j <= z; j++){
							numCRightBound.put(loci[j], numC);
							numCtRightBound.put(loci[j], numCt);
						}
					}
					
				//	System.err.println("right: " + slidingWindowRight.getLength() + "\t" + data.position);
				}
				else{
					int numC = slidingWindowRight.getCReadsNum();
					int numCt = slidingWindowRight.getCtReadsNum();
					LinkedList<Datapoint> tmpValue = slidingWindowRight.addLast(data, true);
					numCLeftBound.put(data.position, slidingWindowRight.getCReadsNum());
					numCtLeftBound.put(data.position, slidingWindowRight.getCtReadsNum());
					numCRightBound.put(data.position, 0);
					numCtRightBound.put(data.position, 0);
					for(Datapoint tmpData: tmpValue){
						
						numCLeftBound.put(tmpData.position, slidingWindowLeft.getCReadsNum());
						numCtLeftBound.put(tmpData.position, slidingWindowLeft.getCtReadsNum());
						slidingWindowLeft.addLast(tmpData, true);
						numCRightBound.put(tmpData.position, numC);
						numCtRightBound.put(tmpData.position, numCt);
						//numC -= tmpData.numC;
						//numCt -= tmpData.numCT;
						

						
					}
					// first add CT reads number when new value add into slidingWindowRight, it will be updated when they are popped out.
					
					
				//	System.err.println("left: " + slidingWindowLeft.getLength() + "\t" + data.position + "\t" + slidingWindowLeft.getCReadsNum() + "\t" + slidingWindowLeft.getCtReadsNum() + "\t" + numC);
				//	System.err.println("right: " + slidingWindowRight.getLength() + "\t" + data.position + "\t" + slidingWindowRight.getCReadsNum() + "\t" + slidingWindowRight.getCtReadsNum() + "\t" + numCt);
				}
			}
			for(GenomeLoc loc : numCtRightBound.keySet()){
				Integer ct = numCtRightBound.get(loc);
				Integer c = numCRightBound.get(loc);
				if(ct < c || c <0 || ct <0)
					System.err.println(c + "\t" + ct + "\t" + loc);
			}
			
			String chr = null;
			int start = -1;
			int end = -1;
			double score = 0;
			int dataPoint = 0; //number of GCH in the segments
			int numC_npr = 0;
			int numT_npr = 0;
			//double sumMethy_Gch_back = 0;
			//int numGch_back = 0;
			//int numC_back = 0;
			//int numT_back = 0;
			
			GenomeLoc startLoc = null;
			GenomeLoc preLoc = null;
			int preState = -1; //record the previous loci's state
			
			double preWeight = -1;
			double postWeight = -1;
			
			for(int i=0; i < hiddenState.length; i++){
				if(preState == -1){
					preState = hiddenState[i]; 
					if(preState == nprState){
						chr = loci[i].getContig();
						//start = loci[i].getStart()-1;
						start = loci[i].getStart()-1;
						startLoc = loci[i];
						score = methyState[i];
						dataPoint=1;
						numC_npr = numCState[i];
						numT_npr = numCTState[i] - numCState[i];

					}	
					
				}
				else{
					if(preState != nprState){
						if(hiddenState[i] == nprState){
								chr = loci[i].getContig();
								preWeight = 1 - nfbsc.getAlpha(hmm, i, preState, hiddenState[i]);
								//double weight = hmm.getAij(preState, hiddenState[i]);
								//start =(int)((loci[i].getStart() + preLoc.getStart())/2) ;
								start =(int)( (loci[i].getStart() - preLoc.getStart()) * preWeight + preLoc.getStart()) ;
								//start =preLoc.getStart() ;
								startLoc = loci[i];
								numC_npr = numCState[i];
								numT_npr = numCTState[i] - numCState[i];
								score = methyState[i];
								dataPoint = 1;
						}

					}
					else{
						if(hiddenState[i] == nprState){
							numC_npr += numCState[i];
							numT_npr += (numCTState[i] - numCState[i]);
							score += methyState[i];
							dataPoint++;
						}
						else{
						//	System.err.println(preLoc);
							//end = (int)((loci[i].getStart() + preLoc.getStart())/2);
							//if(i < hiddenState.length -1){
								postWeight = 1 - nfbsc.getAlpha(hmm, i, preState, hiddenState[i]);
							//}
							//else{
							//	postWeight = 1;
							//}
							//double weight = hmm.getAij(preState, hiddenState[i]);
							end =(int)( (loci[i].getStart() - preLoc.getStart()) * postWeight + preLoc.getStart()) ;
							//end =preLoc.getStart();
							//int numGch_back = backGround.numCLeftBound.size();
							
						//	System.err.println(startLoc);
						//	System.err.println(numCLeftBound.size());
						//	System.err.println(numCRightBound.size());
						//	System.err.println(numCLeftBound.get(startLoc));
						//	System.err.println(numCtLeftBound.get(startLoc));
						//	System.err.println(numCRightBound.get(preLoc));
						//	System.err.println(numCtRightBound.get(preLoc)); 
							int numC_back = numCLeftBound.get(startLoc) + numCRightBound.get(preLoc);
							int numCT_back = numCtLeftBound.get(startLoc) + numCtRightBound.get(preLoc);
							double pValue = getBinomialSigTest(numC_npr, (numC_npr + numT_npr), (double)numC_back/(double)numCT_back, reverseP);
							List<Object> tmp = new LinkedList<Object>();
							
							tmp.add(score);
							tmp.add(dataPoint);
							tmp.add(numC_npr);
							tmp.add(numT_npr);
							tmp.add(numC_back);
							tmp.add(numCT_back);
							//tmp.add(sumMethy_Gch_back);
							//tmp.add(numGch_back);
							tmp.add(end-start);
							tmp.add(1-preWeight);
							tmp.add(1-postWeight);
							tmp.add(pValue);
						//	if((numC_npr + numT_npr >= minCT) && (numC_back + numT_back) >= minCT){
							if((numC_npr + numT_npr) >= minCT && numCT_back >= minCT){
								bedObject bedLine = new bedObject(chr, start, end, Strand.NONE, (List)tmp);
								writer.add(bedLine);
							}
						//	}
						}
					}
					preState = hiddenState[i];
				}
				preLoc = loci[i];
			}
				
					
		}
		
		
		private void getNPRSegmentByHalfLocalCompar(int[] hiddenState, double[] methyState, GenomeLoc[] loci, int[] numCTState, int[] numCState, int nprState, bedObjectWriterImp writer, boolean reverseP, Hmm<ObservationReal> hmm){
			//hash the ct reads in the adjacent window information
			
		
			HashMap<GenomeLoc, Integer> numCtLeftBound = new HashMap<GenomeLoc, Integer>();
			HashMap<GenomeLoc, Integer> numCLeftBound = new HashMap<GenomeLoc, Integer>();
			HashMap<GenomeLoc, Integer> numCtRightBound = new HashMap<GenomeLoc, Integer>();
			HashMap<GenomeLoc, Integer> numCRightBound = new HashMap<GenomeLoc, Integer>();
			//System.err.println(numCTState.length + "\t" + numCState.length);
			SlidingWindow slidingWindowLeft = new SlidingWindow(ctInWindow, gchInWindow, window, nprState);
			SlidingWindow slidingWindowRight = new SlidingWindow(ctInWindow, gchInWindow, window, nprState);
			
			for(int z=0; z < hiddenState.length; z++){
				//System.err.println(z + "\t" + loci.length + "\t" + methyState.length + "\t" + numCTState.length + "\t" + numCState.length + "\t" + hiddenState.length);
				Datapoint data = new Datapoint(loci[z],new ObservationReal(methyState[z]),numCTState[z],numCState[z],hiddenState[z]);
			///	if(loci[z].getStart() == 2007766)
			//		System.err.println(numCRightBound.size() + "\t" + slidingWindowRight.getCReadsNum() + "\t" + slidingWindowRight.getLength() + "\t" + slidingWindowLeft.getLength() + "\t" + z);
				if(slidingWindowLeft.getLength() < window || slidingWindowLeft.getGchNum() < gchInWindow || slidingWindowLeft.getCtReadsNum() < ctInWindow){ // in the beginning of the chromosome
					slidingWindowLeft.addLast(data);
					if(z==hiddenState.length-1){
						int numC = slidingWindowLeft.getCReadsNum();
						int numCt = slidingWindowLeft.getCtReadsNum();
						for(int j = 0; j <= z; j++){
							numCLeftBound.put(loci[j], numC);
							numCtLeftBound.put(loci[j], numCt);
							numCRightBound.put(loci[j], 0);
							numCtRightBound.put(loci[j], 0);
						}
					}
					
				//	System.err.println("left: " + slidingWindowLeft.getLength() + "\t" + data.position);
				}
				else if(slidingWindowRight.windowList.isEmpty()){
					int numC = slidingWindowLeft.getCReadsNum();
					int numCt = slidingWindowLeft.getCtReadsNum();
					for(int j = 0; j < z; j++){
						numCLeftBound.put(loci[j], numC);
						numCtLeftBound.put(loci[j], numCt);
						numCRightBound.put(loci[j], 0);
						numCtRightBound.put(loci[j], 0);
					}
					slidingWindowRight.addLast(data);
					numCLeftBound.put(data.position, numC);
					numCtLeftBound.put(data.position, numCt);
					numCRightBound.put(data.position, 0);
					numCtRightBound.put(data.position, 0);
				}
				else if(slidingWindowRight.getLength() < window || slidingWindowRight.getGchNum() < gchInWindow || slidingWindowRight.getCtReadsNum() < ctInWindow){
					slidingWindowRight.addLast(data);
					numCLeftBound.put(data.position, slidingWindowRight.getCReadsNum());
					numCtLeftBound.put(data.position, slidingWindowRight.getCtReadsNum());
					numCRightBound.put(data.position, 0);
					numCtRightBound.put(data.position, 0);
					if(z==hiddenState.length-1){
						int numC = slidingWindowRight.getCReadsNum();
						int numCt = slidingWindowRight.getCtReadsNum();
						for(int j = z - slidingWindowRight.windowList.size(); j <= z; j++){
							numCRightBound.put(loci[j], numC);
							numCtRightBound.put(loci[j], numCt);
						}
					}
					
				//	System.err.println("right: " + slidingWindowRight.getLength() + "\t" + data.position);
				}
				else{
					int numC = slidingWindowRight.getCReadsNum();
					int numCt = slidingWindowRight.getCtReadsNum();
					LinkedList<Datapoint> tmpValue = slidingWindowRight.addLast(data, true);
					numCLeftBound.put(data.position, slidingWindowRight.getCReadsNum());
					numCtLeftBound.put(data.position, slidingWindowRight.getCtReadsNum());
					numCRightBound.put(data.position, 0);
					numCtRightBound.put(data.position, 0);
					for(Datapoint tmpData: tmpValue){
						
						numCLeftBound.put(tmpData.position, slidingWindowLeft.getCReadsNum());
						numCtLeftBound.put(tmpData.position, slidingWindowLeft.getCtReadsNum());
						slidingWindowLeft.addLast(tmpData, true);
						numCRightBound.put(tmpData.position, numC);
						numCtRightBound.put(tmpData.position, numCt);
						//numC -= tmpData.numC;
						//numCt -= tmpData.numCT;
						

						
					}
					// first add CT reads number when new value add into slidingWindowRight, it will be updated when they are popped out.
					
					
				//	System.err.println("left: " + slidingWindowLeft.getLength() + "\t" + data.position + "\t" + slidingWindowLeft.getCReadsNum() + "\t" + slidingWindowLeft.getCtReadsNum() + "\t" + numC);
				//	System.err.println("right: " + slidingWindowRight.getLength() + "\t" + data.position + "\t" + slidingWindowRight.getCReadsNum() + "\t" + slidingWindowRight.getCtReadsNum() + "\t" + numCt);
				}
			}
			for(GenomeLoc loc : numCtRightBound.keySet()){
				Integer ct = numCtRightBound.get(loc);
				Integer c = numCRightBound.get(loc);
				if(ct < c || c <0 || ct <0)
					System.err.println(c + "\t" + ct + "\t" + loc);
			}
			
			String chr = null;
			int start = -1;
			int end = -1;
			double score = 0;
			int dataPoint = 0; //number of GCH in the segments
			int numC_npr = 0;
			int numT_npr = 0;
			//double sumMethy_Gch_back = 0;
			//int numGch_back = 0;
			//int numC_back = 0;
			//int numT_back = 0;
			
			GenomeLoc startLoc = null;
			GenomeLoc preLoc = null;
			int preState = -1; //record the previous loci's state
			
			double preWeight = -1;
			double postWeight = -1;
			
			for(int i=0; i < hiddenState.length; i++){
				if(preState == -1){
					preState = hiddenState[i]; 
					if(preState == nprState){
						chr = loci[i].getContig();
						//start = loci[i].getStart()-1;
						start = loci[i].getStart()-1;
						startLoc = loci[i];
						score = methyState[i];
						dataPoint=1;
						numC_npr = numCState[i];
						numT_npr = numCTState[i] - numCState[i];

					}	
					
				}
				else{
					if(preState != nprState){
						if(hiddenState[i] == nprState){
								chr = loci[i].getContig();
								double weight = hmm.getAij(preState, hiddenState[i]);
								//start =(int)((loci[i].getStart() + preLoc.getStart())/2) ;
								start =(int)( (loci[i].getStart() - preLoc.getStart()) * weight + preLoc.getStart()) ;
								startLoc = loci[i];
								numC_npr = numCState[i];
								numT_npr = numCTState[i] - numCState[i];
								score = methyState[i];
								dataPoint = 1;
						}

					}
					else{
						if(hiddenState[i] == nprState){
							numC_npr += numCState[i];
							numT_npr += (numCTState[i] - numCState[i]);
							score += methyState[i];
							dataPoint++;
						}
						else{
						//	System.err.println(preLoc);
							//end = (int)((loci[i].getStart() + preLoc.getStart())/2);
							double weight = hmm.getAij(preState, hiddenState[i]);
							end =(int)( (loci[i].getStart() - preLoc.getStart()) * weight + preLoc.getStart()) ;
							//int numGch_back = backGround.numCLeftBound.size();
							
	
							int numC_back = numCLeftBound.get(startLoc) + numCRightBound.get(preLoc);
							int numCT_back = numCtLeftBound.get(startLoc) + numCtRightBound.get(preLoc);
							double pValue = getBinomialSigTest(numC_npr, (numC_npr + numT_npr), (double)numC_back/(double)numCT_back, reverseP);
							List<Object> tmp = new LinkedList<Object>();
							
							tmp.add(score);
							tmp.add(dataPoint);
							tmp.add(numC_npr);
							tmp.add(numT_npr);
							tmp.add(numC_back);
							tmp.add(numCT_back);
							//tmp.add(sumMethy_Gch_back);
							//tmp.add(numGch_back);
							tmp.add(end-start);
							tmp.add(pValue);
						//	if((numC_npr + numT_npr >= minCT) && (numC_back + numT_back) >= minCT){
							if((numC_npr + numT_npr) >= minCT && numCT_back >= minCT){
								bedObject bedLine = new bedObject(chr, start, end, Strand.NONE, (List)tmp);
								writer.add(bedLine);
							}
						//	}
						}
					}
					preState = hiddenState[i];
				}
				preLoc = loci[i];
			}
				
					
		}
		

	private double getBinomialSigTest(int k, int n, double pSucess, boolean reverseP){
		
		//double p;
		//if(k==0){
		//	p = Probability.binomial(k,n,pValue);
		//}
		//else{
		//	p = Probability.binomial(k,n,pValue)-Probability.binomial(k-1,n,pValue);
		//}
		//System.err.println(k + "\t" + n + "\t" + pValue + "\t" + p);
		BinomialDistributionImpl binomial = new BinomialDistributionImpl(n, pSucess);
		
		//return binomial.probability(k);
		 double p = Double.NaN;
		try {
			if(reverseP){
				p = 1-binomial.cumulativeProbability(k)+binomial.probability(k); //p(X>= x)
			}
			else{
				p = binomial.cumulativeProbability(k); //p(X<= x)
			}
			
		} catch (MathException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return p;
		//
		
	}
	
	private double getFisherPvalue(int methy, int unmethy, int methyAdj, int unmethyAdj, boolean rightTail){
		//	System.err.println(methy + "\t" + unmethy + "\t" + methyAdj + "\t" + unmethyAdj);
		//	ContingencyTable2x2 table = new ContingencyTable2x2(methy, methyAdj, unmethy, unmethyAdj);
		//	ContingencyTable2x2 table = new ContingencyTable2x2(22, 0, expectBaseCount, expectOtherBaseCount);
		//	FishersExactTest fisherExact = new FishersExactTest(table, alternative);
		//	return fisherExact.getSP();
		double pValue = Double.NaN;
		//System.err.println(methy + "\t" + unmethy + "\t" + methyAdj + "\t" + unmethyAdj);
		if(methy+unmethy == 0 || unmethyAdj + methyAdj == 0){
			return pValue;
		}
		//if(methy + unmethy + methyAdj + unmethyAdj >= 10000){
		//	System.err.println(methy + "\t" + unmethy + "\t" + methyAdj + "\t" + unmethyAdj);
		//	ContingencyTable2x2 table = new ContingencyTable2x2(methy,methyAdj, unmethy,  unmethyAdj);
			
		//	ChiSquaredTest chiSquaredTest = new ChiSquaredTest(table);
		//	pValue = chiSquaredTest.getSP();
		//}
		//else{
		int sum = methy + methyAdj + unmethy + unmethyAdj + 100;
			FisherExactTest fisherExact = new FisherExactTest(sum);
			if(rightTail){
				pValue = fisherExact.getRightTailedP(methy,methyAdj, unmethy,  unmethyAdj);
			}
			else{
				pValue = fisherExact.getLeftTailedP(methy,methyAdj, unmethy,  unmethyAdj);
			}
			 
	//	}
			
			// System.err.println(22 + "\t" + 0 + "\t" + expectBaseCount + "\t" + expectOtherBaseCount + "\t" + pValue);
			  return pValue;
		}
	
	public class Datum{		
		public ArrayList<LinkedList<ObservationReal>> value;
		public ArrayList<LinkedList<GenomeLoc>> position;
		public ArrayList<LinkedList<Integer>> numCT;
		public ArrayList<LinkedList<Integer>> numC;
		public ArrayList<LinkedList<ObservationMethy>> methy;
//		public HashMap<GenomeLoc, Integer> numCtLeftBound;
//		public HashMap<GenomeLoc, Integer> numCLeftBound;
//		public HashMap<GenomeLoc, Integer> numCtRightBound;
//		public HashMap<GenomeLoc, Integer> numCRightBound;
	
	}
	
	public class Datapoint{
		public ObservationReal value;
		public GenomeLoc position;
		public ObservationMethy methy;
		public int numCT;
		public int numC;
		public int hmmState;

		public Datapoint(GenomeLoc position, ObservationReal value, int numCT){
			this.position = position;
			this.value = value;
			this.numCT = numCT;
			
		}
		
		public Datapoint(GenomeLoc position, ObservationReal value, int numCT, int numC){
			this.position = position;
			this.value = value;
			this.numCT = numCT;
			this.numC = numC;
			this.methy = new ObservationMethy(value.value);
			methy.setCoverage(numCT);
		}
		
		public Datapoint(GenomeLoc position, ObservationReal value, int numCT, int numC, int hmmState){
			this.position = position;
			this.value = value;
			this.numCT = numCT;
			this.numC = numC;
			this.hmmState = hmmState;
		}
		
	}
	

	
	public class SlidingWindow{
		public LinkedList<Datapoint> windowList;
		private int numC = 0;
		private int numCT = 0;
		private int minCT;
		private int minGch;
		private long minLen;
		private int segHmmState;
		
		
		public SlidingWindow(int minCT, int minGch, long minLen, int segHmmState){
			windowList = new LinkedList<Datapoint>();
			this.minCT = minCT;
			this.minGch = minGch;
			this.minLen = minLen;
			this.segHmmState = segHmmState;
			
		}
		
		public void addLast(Datapoint data){
			windowList.offerLast(data);
			if(data.hmmState != segHmmState){
				numCT += data.numCT;
				numC += data.numC;
			}
			
			
		}
		
		public void addFirst(Datapoint data){
			windowList.offerFirst(data);
			if(data.hmmState != segHmmState){
				numCT += data.numCT;
				numC += data.numC;
			}
			
			
		}
		
		public LinkedList<Datapoint> addLast(Datapoint data, boolean automateRemoveFirstBatch){
			LinkedList<Datapoint> dataList = new LinkedList<Datapoint>();
			addLast(data);
			while(getGchNum() > minGch && getLength() > minLen && getCtReadsNum() > minCT){
				dataList.offerLast(removeFirst());
			}
			if(!dataList.isEmpty()){
				addFirst(dataList.pollLast());
			}
			
			return dataList;
		}
		
		public Datapoint removeLast(){
			Datapoint data = windowList.pollLast();
			if(data.hmmState != segHmmState){
				numCT -= data.numCT;
				numC -= data.numC;
			}
			
			return data;
		}
		
		public Datapoint removeFirst(){
			Datapoint data = windowList.pollFirst();
			if(data.hmmState != segHmmState){
				numCT -= data.numCT;
				numC -= data.numC;
			}
			
			return data;
		}
		
		public double getMean(){
			
			return (double)numC/(double)numCT;
		}
		
		public int getGchNum(){
			return windowList.size();
		}
		
		public int getLength(){
			if(!windowList.isEmpty()){
				return windowList.peekLast().position.distance(windowList.peekFirst().position);
			}
			return 0;
		}
		
		public int getCtReadsNum(){

			return numCT;
		}
		
		public int getCReadsNum(){
			
			return numC;
		}
		
	}
	
}
