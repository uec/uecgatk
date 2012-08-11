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
import java.util.LinkedList;
import java.util.List;
import java.util.Iterator;


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
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import edu.usc.epigenome.uecgatk.YapingWriter.bedObject;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObjectWriterImp;
import edu.usc.epigenome.uecgatk.distribution.OpdfBeta;
import edu.usc.epigenome.uecgatk.distribution.OpdfBetaBinomialFactory;
import edu.usc.epigenome.uecgatk.distribution.OpdfBetaFactory;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteVCFConstants;
import edu.usc.epigenome.uecgatk.NOMeSeqWalker.OpdfBetaWriter;
import edu.usc.epigenome.uecgatk.NOMeSeqWalker.OpdfBetaReader;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationDiscrete;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.OpdfDiscrete;
import be.ac.ulg.montefiore.run.jahmm.OpdfDiscreteFactory;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianFactory;
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
	
	@Output(fullName = "MPR_segment_file", shortName = "mprSeg", doc = "write segment of MPR/NPR region into a file (no training mode)", required = false)
	public String mprSegFile = null;
	
	@Argument(fullName = "min_ct_coverage", shortName = "minCT", doc = "minimum number of CT reads for count methylation level, default: 1", required = false)
	public int minCT = 1;
	
	@Argument(fullName = "training_mode", shortName = "train", doc = "enable the training mode and output HMM parameter, default: not enabled", required = false)
	public boolean train = false;
	
	@Argument(fullName = "hmm_file", shortName = "hmm", doc = "read/write HMM model from/to a file, default: read HMM parameters from a file", required = true)
	public String hmmFile = null;
	
	@Argument(fullName = "max_gap_size", shortName = "gap", doc = "max gap size, default: 1000", required = false)
	public int gap = 1000;
	
	@Argument(fullName = "min_data_point", shortName = "dataP", doc = "minimum data point, default: 5", required = false)
	public int dataP = 2;
	
	@Argument(fullName = "num_of_states", shortName = "states", doc = "number of states, default: 2", required = false)
	public static int states = 2;
	
	private static int STATES=states;
	
	private double TOLERENCE=1e-5;
	
	private int MAXIMUM_GAP_SIZE=gap; // maximum gap allowed for split different sequence for training & decoding
	
	private int MAXIMUM_DATA_POINTS=1000000000;
	
	private int MINIMUM_DATA_POINTS=dataP;
	
	private Hmm<ObservationReal> hmm = null;
	
	private bedObjectWriterImp bedWriter = null;
	
	private bedObjectWriterImp segWriter = null;
	
	private bedObjectWriterImp nprSegWriter = null;
	
	
	
	public void initialize(){
		if(!train){
			bedWriter = new bedObjectWriterImp(new File(resultFile));
			segWriter = new bedObjectWriterImp(new File(segmentFile));
			nprSegWriter = new bedObjectWriterImp(new File(mprSegFile));
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
								methyValue += (Math.random())/100000;
							else if(numT == 0)
								methyValue -= (Math.random())/100000;
							//double methyValue =(double) numC / (double) (numC + numT);
							//System.err.println(ref.getLocus());
							//System.err.println(new ObservationReal(methyValue));
							Datapoint dat = new Datapoint(ref.getLocus(),new ObservationReal(methyValue), numC + numT);
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
		tmp.position = new ArrayList<LinkedList<GenomeLoc>>();
		tmp.value = new ArrayList<LinkedList<ObservationReal>>();
		tmp.numCT = new ArrayList<LinkedList<Integer>>();
		tmp.position.add(position);
		tmp.value.add(value);
		tmp.numCT.add(numCT);
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
		LinkedList<Integer> numCT = sum.numCT.remove(sum.numCT.size()-1);
		if(!position.isEmpty()){
			if( position.size() >= MAXIMUM_DATA_POINTS|| position.peekLast().distance(value.position) >= MAXIMUM_GAP_SIZE || !position.peekLast().onSameContig(value.position)){
				if(position.size() >= MINIMUM_DATA_POINTS){
					sum.position.add(position);
					sum.value.add(data);
					sum.numCT.add(numCT);
				}
				
				LinkedList<GenomeLoc> newPosition = new LinkedList<GenomeLoc>();
				LinkedList<ObservationReal> newValue = new LinkedList<ObservationReal>();
				LinkedList<Integer> newNumCT = new LinkedList<Integer>();
				newPosition.offerLast(value.position);
				newValue.offerLast(value.value);
				newNumCT.offerLast(value.numCT);
				sum.position.add(newPosition);
				sum.value.add(newValue);
				sum.numCT.add(newNumCT);
				return sum;
			}
		}
		position.offerLast(value.position);
		data.offerLast(value.value);
		numCT.offerLast(value.numCT);
		sum.position.add(position);
		sum.value.add(data);
		sum.numCT.add(numCT);
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
	//	ArrayList<List<ObservationReal>> values = new ArrayList<List<ObservationReal>>();
		System.out.println("sequence data size: " + result.value.size());
		int z = 0;
		for(LinkedList<ObservationReal> tmp : result.value){
			System.out.println("size: " + tmp.size());
			if(tmp.size() < dataP){
				result.value.remove(z);
				result.position.remove(z);
				result.numCT.remove(z);
			}	
			z++;
		}
	//	values.add(result.value);
		if(train){
			System.out.println("training....");
			hmm = buildInitHmmByBeta(result.value);
			BaumWelchScaledLearner bwl = new BaumWelchScaledLearner();

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
			while(Math.abs(distance) >= TOLERENCE){
				i++;
				
				try {
					prevHmm = hmm.clone();
				} catch (CloneNotSupportedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				System.out.println("HMM pre:\n" + prevHmm);
				
				hmm = bwl.iterate(hmm, result.value);
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
			int ndrState = 0;
			double maxMean = 0;
			int nprState = 1;
			double minMean = 1.1;
			//suppose NDR region have the highest mean GCH methyation value, NPR region have lowest mean GCH methylation value
			for(int m = 0; m < hmm.nbStates(); m++){
				if(((OpdfBeta)hmm.getOpdf(m)).mean() > maxMean){
					maxMean = ((OpdfBeta)hmm.getOpdf(m)).mean();
					ndrState = m;
				}
				if(((OpdfBeta)hmm.getOpdf(m)).mean() < minMean){
					minMean = ((OpdfBeta)hmm.getOpdf(m)).mean();
					nprState = m;
				}
			}
			
			double[] randomSeqs = getPermutatedSeqs(result.value, hmm, ndrState);
			
			int j=0;
			
			for(LinkedList<ObservationReal> value : result.value){
				int[] hiddenState = hmm.mostLikelyStateSequence(value);
				
				GenomeLoc[] loci = new GenomeLoc[result.position.get(j).size()];
				Iterator<GenomeLoc> it = result.position.get(j).iterator();
				int ii=0;
				while(it.hasNext()){
					loci[ii] = it.next();
					ii++;
				}
				
				ii=0;
				double[] methyState = new double[value.size()];
				Iterator<ObservationReal> it2 = value.iterator();
				while(it2.hasNext()){
					methyState[ii] = it2.next().value;
					ii++;
				}
				
				ii=0;
				int[] numCTState = new int[result.numCT.get(j).size()];
				Iterator<Integer> it3 = result.numCT.get(j).iterator();
				while(it3.hasNext()){
					numCTState[ii] = it3.next();
					ii++;
				}
				//getNDRSegment(hiddenState, methyState, loci, numCTState, randomSeqs, ndrState);
				getNPRSegmentByLocalCompar(hiddenState, methyState, loci, numCTState, ndrState, segWriter);
				getNPRSegmentByLocalCompar(hiddenState, methyState, loci, numCTState, nprState, nprSegWriter);
				for(int i = 0; i < hiddenState.length; i++){
					List<Object> tmp = new LinkedList<Object>();
					tmp.add(hiddenState[i]);
					tmp.add(methyState[i]);
					tmp.add(numCTState[i]);
					bedObject bedLine = new bedObject(loci[i].getContig(), loci[i].getStart()-1, loci[i].getStop(), Strand.NONE, (List)tmp);
					bedWriter.add(bedLine);
				}
				j++;
			}

			bedWriter.close();
			segWriter.close();
			nprSegWriter.close();
		}
	}
	
	private static Hmm<ObservationReal> buildInitHmm(ArrayList<LinkedList<ObservationReal>> seqs)
	{	

		KMeansLearner<ObservationReal> kl = new KMeansLearner<ObservationReal>(STATES, new OpdfGaussianFactory(),
				seqs);
		System.out.println("KMeansLearner...");
		Hmm<ObservationReal> hmm = kl.learn();
		return hmm;
	}
	
	private static Hmm<ObservationReal> buildInitHmmByBeta(ArrayList<LinkedList<ObservationReal>> seqs)
	{	

		KMeansLearner<ObservationReal> kl = new KMeansLearner<ObservationReal>(STATES, new OpdfBetaFactory(),
				seqs);
		System.out.println("KMeansLearner...");
		Hmm<ObservationReal> hmm = kl.learn();
		return hmm;
	}
	
	private static Hmm<ObservationReal> buildInitHmmRandomlyByBeta(ArrayList<LinkedList<ObservationReal>> seqs)
	{	

	//	KMeansLearner<ObservationReal> kl = new KMeansLearner<ObservationReal>(STATES, new OpdfBetaFactory(),
		//		seqs);
	//	System.out.println("KMeansLearner...");
	//	Hmm<ObservationReal> hmm = kl.learn();
		
		/*
		Hmm<ObservationReal> hmm = 
				new Hmm<ObservationReal>(2,new OpdfBetaFactory());
		double tmp3 = Math.random();
			hmm.setPi(0, tmp3);
			hmm.setPi(1, 1-tmp3);
			
			hmm.setOpdf(0, new OpdfBeta(Math.random(),Math.random()));
			hmm.setOpdf(1, new OpdfBeta(Math.random(),Math.random()));
			double tmp1 = Math.random();
			
			hmm.setAij(0, 1, tmp1);
			hmm.setAij(0, 0, 1-tmp1);
			double tmp2 = Math.random();
			hmm.setAij(1, 0, tmp1);
			hmm.setAij(1, 1, 1-tmp1);
		*/
		Hmm<ObservationReal> hmm = 
				new Hmm<ObservationReal>(2,new OpdfBetaFactory());
		double tmp3 = Math.random();
			hmm.setPi(0, 0.7);
			hmm.setPi(1, 0.3);
			
			hmm.setOpdf(0, new OpdfBeta(0.3,1.1));
			hmm.setOpdf(1, new OpdfBeta(1.0,0.3));
			double tmp1 = Math.random();
			
			hmm.setAij(0, 1, 0.2);
			hmm.setAij(0, 0, 0.8);
			double tmp2 = Math.random();
			hmm.setAij(1, 0, 0.1);
			hmm.setAij(1, 1, 0.9);
		return hmm;
	}
	
	
	private static Hmm<ObservationMethy> buildInitHmmByBetaBinomial(ArrayList<LinkedList<ObservationMethy>> seqs)
	{	

		KMeansLearner<ObservationMethy> kl = new KMeansLearner<ObservationMethy>(STATES, new OpdfBetaBinomialFactory(),
				seqs);
		System.out.println("KMeansLearner...");
		Hmm<ObservationMethy> hmm = kl.learn();
		return hmm;
	}
	

	private double[] getPermutatedSeqs(ArrayList<LinkedList<ObservationReal>> seqs, Hmm<ObservationReal> hmm, int ndrState){
		ArrayList<ObservationReal> randomValueList = new ArrayList<ObservationReal>();
		
		Iterator<LinkedList<ObservationReal>> it = seqs.iterator();
		while(it.hasNext()){
			randomValueList.addAll(it.next());
		}
		Collections.shuffle(randomValueList);
		int[] hiddenState = hmm.mostLikelyStateSequence(randomValueList);
		double[] score = getNDRSegment(hiddenState, randomValueList, ndrState);
		return score;
	}
	
	private double[] getNDRSegment(int[] hiddenState, ArrayList<ObservationReal> randomValueList, int ndrState){
		ArrayList<Double> scoreList = new ArrayList<Double>();
		double score = 0;
		int preState = -1;
		int i = 0;
		Iterator<ObservationReal> it = randomValueList.iterator();
		while(it.hasNext()){
			double methyState = it.next().value;
			if(preState == -1){
				preState = hiddenState[i];
				if(preState == ndrState){
					
					score += methyState;
				}
				i++;
				continue;
			}
			if(preState != ndrState && hiddenState[i] == ndrState){

				score += methyState;
				preState = hiddenState[i];
			}
			else if(preState == ndrState){
				if(hiddenState[i] != ndrState){
					preState = hiddenState[i];
					scoreList.add(score);
					score = 0;
				}
				else{
					preState = hiddenState[i];
					score += methyState;
				}
				
			}
			i++;
		}
		double[] scoreCollection = new double[scoreList.size()];
		Iterator<Double> it2 = scoreList.iterator();
		i = 0;
		while(it2.hasNext()){
			scoreCollection[i++] = it2.next();
		}
		return scoreCollection;
	}
	
	//NDR segement: chr, start, end, strand, score(methylation sum of this segment), GCH_number, numCT (total number of CT in this segment), seg_length, raw_p_value 
	private void getNDRSegment(int[] hiddenState, double[] methyState, GenomeLoc[] loci, int[] numCTState, double[] randomSeqs, int ndrState){
		String chr = null;
		int start = -1;
		int end = -1;
		double score = 0;
		int preState = -1;
		GenomeLoc preLoc = null;
		int dataPoint = 0;
		int numCT = 0;
		for(int i=0; i < hiddenState.length; i++){
			if(preState == -1){
				preState = hiddenState[i];
				preLoc = loci[i]; 
				if(preState == ndrState){
					chr = loci[i].getContig();
					start = loci[i].getStart()-1;
					score += methyState[i];
					dataPoint++;
					numCT += numCTState[i];
				}	
				continue;
			}
			if(preState != ndrState && hiddenState[i] == ndrState){
				chr = loci[i].getContig();
				start = loci[i].getStart()-1;
				score += methyState[i];
				numCT += numCTState[i];
				dataPoint++;
				preState = hiddenState[i];
				preLoc = loci[i];
			}
			else if(preState == ndrState){
				if(hiddenState[i] != ndrState){
					end = preLoc.getStart();
					preState = hiddenState[i];
					List<Object> tmp = new LinkedList<Object>();
					tmp.add(score);
					tmp.add(dataPoint);
					tmp.add(numCT);
					tmp.add(end-start);
					
					tmp.add(getPvalue(score, randomSeqs));
					bedObject bedLine = new bedObject(chr, start, end, Strand.NONE, (List)tmp);
					segWriter.add(bedLine);
					score = 0;
					numCT = 0;
					dataPoint = 0;
					preLoc = loci[i];
				}
				else{
					preState = hiddenState[i];
					preLoc = loci[i];
					score += methyState[i];
					numCT += numCTState[i];
					dataPoint++;
				}
				
			}
			
		}
	}
	
	//output NPR segment format: chr, start, end, strand, score(methylation sum of this segment), GCH_number, numC (total number of C in this segment),
	//numT (total number of T in this segment), numC_adj (total number of C in the adjacent segment), numT_adj (total number of T in the adjacent segment),seg_length, raw_p_value
		private void getNPRSegmentByLocalCompar(int[] hiddenState, double[] methyState, GenomeLoc[] loci, int[] numCTState, int nprState, bedObjectWriterImp writer){
			String chr = null;
			int start = -1;
			int end = -1;
			double score = 0;
			int preState = -1; //record the previous loci's state
			boolean passNprSeg = false;
			int dataPoint = 0; //number of GCH in the segments

			int numC_pre = 0;
			int numT_pre = 0;
			int numC_after = 0;
			int numT_after = 0;
			int numC_npr = 0;
			int numT_npr = 0;
			
			for(int i=0; i < hiddenState.length; i++){
				if(preState == -1){
					preState = hiddenState[i]; 
					if(preState == nprState){
						chr = loci[i].getContig();
						start = loci[i].getStart()-1;
						score = methyState[i];
						dataPoint=1;
						numC_npr = (int)(numCTState[i] * methyState[i]);
						numT_npr = (int)(numCTState[i] * (1-methyState[i]));

					}	
					else{
						numC_pre = (int)(numCTState[i] * methyState[i]);
						numT_pre = (int)(numCTState[i] * (1-methyState[i]));
					}
					continue;
				}
				else{
					if(preState != nprState){
						if(hiddenState[i] == nprState){
							if(passNprSeg){
								double pSuccess = (double)(numC_pre + numT_pre)/(double)(numC_pre + numT_pre + numC_after + numT_after);
								double pValue = getBinomialSigTest(numC_npr, (numC_npr+numT_npr), pSuccess);
								List<Object> tmp = new LinkedList<Object>();
								tmp.add(score);
								tmp.add(dataPoint);
								tmp.add(numC_npr);
								tmp.add(numT_npr);
								tmp.add((numC_pre + numC_after));
								tmp.add((numT_pre + numT_after));
								tmp.add(end-start);
								
								tmp.add(pValue);
								bedObject bedLine = new bedObject(chr, start, end, Strand.NONE, (List)tmp);
								writer.add(bedLine);
								chr = loci[i].getContig();
								start = loci[i].getStart()-1;
								end = loci[i].getStart();
								score = methyState[i];
								dataPoint = 1;
								numC_npr = (int)(numCTState[i] * methyState[i]);
								numT_npr = (int)(numCTState[i] * (1-methyState[i]));
								numC_pre = numC_after;
								numT_pre = numT_after;
								numC_after = 0;
								numT_after = 0;
								
							}
							else{
								chr = loci[i].getContig();
								start = loci[i].getStart()-1;
								numC_npr += (int)(numCTState[i] * methyState[i]);
								numT_npr += (int)(numCTState[i] * (1-methyState[i]));
								score += methyState[i];
								dataPoint++;
							}
						}
						else{
							if(passNprSeg){
								numC_after += (int)(numCTState[i] * methyState[i]);
								numT_after += (int)(numCTState[i] * (1-methyState[i]));
							}
							else{
								numC_pre += (int)(numCTState[i] * methyState[i]);
								numT_pre += (int)(numCTState[i] * (1-methyState[i]));
							}
						}
						
					}
					else{
						if(hiddenState[i] == nprState){
							numC_npr += (int)(numCTState[i] * methyState[i]);
							numT_npr += (int)(numCTState[i] * (1-methyState[i]));
							score += methyState[i];
							dataPoint++;
						}
						else{
							passNprSeg = true;
							end = loci[i].getStart();
							numC_after += (int)(numCTState[i] * methyState[i]);
							numT_after += (int)(numCTState[i] * (1-methyState[i]));
							
						}
					}
				}
			}
		}
	
	private double getPvalue(double value, double[] distribution){
		int rank = 0;
		for(int i = 0; i < distribution.length; i++){
			if(distribution[i] >= value){
				rank++;
			}
		}
		return (double)rank/(double)distribution.length;
	}
	
	private double getBinomialSigTest(int k, int n, double pSucess){
		
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
			p = binomial.cumulativeProbability(k); //p(X<= x)
		} catch (MathException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return p;
		//
		
	}
	
	public class Datum{
		public ArrayList<LinkedList<ObservationReal>> value;
		public ArrayList<LinkedList<GenomeLoc>> position;
		public ArrayList<LinkedList<Integer>> numCT;
	}
	
	public class Datapoint{
		public ObservationReal value;
		public GenomeLoc position;
		public int numCT;
		public Datapoint(GenomeLoc position, ObservationReal value, int numCT){
			this.position = position;
			this.value = value;
			this.numCT = numCT;
		}
	}
	
	
}
