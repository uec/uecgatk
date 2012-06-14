/**
 * 
 */
package edu.usc.epigenome.uecgatk.NOMeSeqWalker;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Downsample;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.GenomeLoc;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.ObservationDiscrete;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.OpdfDiscrete;
import be.ac.ulg.montefiore.run.jahmm.OpdfDiscreteFactory;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianFactory;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianMixture;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianMixtureFactory;
import be.ac.ulg.montefiore.run.jahmm.apps.sample.SimpleExample.Packet;
import be.ac.ulg.montefiore.run.jahmm.draw.GenericHmmDrawerDot;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchLearner;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchScaledLearner;
import be.ac.ulg.montefiore.run.jahmm.learn.KMeansLearner;
import be.ac.ulg.montefiore.run.jahmm.toolbox.KullbackLeiblerDistanceCalculator;
import be.ac.ulg.montefiore.run.jahmm.toolbox.MarkovGenerator;

import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteArgumentCollection;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteGenotyperEngine;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteVariantCallContext;
import edu.usc.epigenome.uecgatk.bisulfiteIndels.BisBAQ;
import edu.usc.epigenome.uecgatk.bisulfiteIndels.BisBAQMode;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jun 5, 2012 8:56:42 PM
 * 
 */

@BisBAQMode(QualityMode = BisBAQ.QualityMode.ADD_TAG, ApplicationTime = BisBAQ.ApplicationTime.ON_INPUT)
@ReadFilters( {UnmappedReadFilter.class, BadMateFilter.class, NotPrimaryAlignmentFilter.class, DuplicateReadFilter.class} ) // Filter out all reads with zero mapping quality
@Reference(window=@Window(start=-500,stop=500))
@By(DataSource.REFERENCE)
@Downsample(by=DownsampleType.NONE)
public class HMMestimateParaWalker extends LocusWalker<BisulfiteVariantCallContext, Long> implements TreeReducible<Long> {

	@ArgumentCollection private static BisulfiteArgumentCollection BAC = new BisulfiteArgumentCollection();
	
	@Argument(fullName = "minium_CT_reads_count", shortName = "minCTdepth", doc = "minium number of CT reads should contained to calculate methylation value", required = false)
    public int minCTdepth = 1;
	
	@Output(doc = "output training process", required = true)
    public PrintStream writer = null;
	
	private LinkedList<ObservationReal> methyStatusListGch = null;
	
	private LinkedList<GenomeLoc> positionListGch = null;
	
	private BisulfiteGenotyperEngine BG_engine = null;
	
	private int emptyLociInList=0;
	
	private int numSeq = 0;
	private int maxLen = 0;
	
	private int MAXIMUM_GAP_SIZE=150;
	
	private int MINIMUM_DATA_POINTS=5;
	
	private static int STATES=2;
	
	private int MAXIMUM_DATA_POINTS=1000000;
	
	private double TOLERENCE=1e-10;
	
	private HiddenMarkov mkv = null;
	
	private ArrayList<List<ObservationReal>> storage;
	
	private ArrayList<GenomeLoc[]> locStorage;
	
	public void initialize(){
		BAC.cytosineContextsAcquired.add("GCH,2");
		 BAC.makeCytosine();
		 methyStatusListGch = new LinkedList<ObservationReal>();
		 positionListGch = new LinkedList<GenomeLoc>();
		 mkv = new HiddenMarkov(2, 10);
		 storage = new ArrayList<List<ObservationReal>>();
		 locStorage = new ArrayList<GenomeLoc[]>();
		// mkv.setNumObSeq(10000000);
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.LocusWalker#map(org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker, org.broadinstitute.sting.gatk.contexts.ReferenceContext, org.broadinstitute.sting.gatk.contexts.AlignmentContext)
	 */
	@Override
	public BisulfiteVariantCallContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

		if(context == null){
			emptyLociInList++;
			if(emptyLociInList >= MAXIMUM_GAP_SIZE || ((!positionListGch.isEmpty()) && positionListGch.peekLast().distance(ref.getLocus()) >= MAXIMUM_GAP_SIZE) || methyStatusListGch.size()>=MAXIMUM_DATA_POINTS){
				//send list to HMM model, get state sequence back
				if(!methyStatusListGch.isEmpty() && methyStatusListGch.size() >= MINIMUM_DATA_POINTS){
					//sendToHmm();
					//System.err.println(ref.getLocus().getContig() + "\t" + ref.getLocus().getStart());
				}
					
				methyStatusListGch.clear();
				positionListGch.clear();
				//System.err.println("gapsize: " + emptyLociInList);
				emptyLociInList=0;
				return null;
			}
		}
		
		BG_engine = new BisulfiteGenotyperEngine(tracker, ref, context, BAC, getToolkit());
 		BisulfiteVariantCallContext bvc = BG_engine.getBisulfiteVariantCallContext();
		
 		if(bvc == null || bvc.getSummaryAcrossRG().cytosinePatternConfirmedSet ==null){
 			emptyLociInList++;
			if(emptyLociInList >= MAXIMUM_GAP_SIZE || ((!positionListGch.isEmpty()) && positionListGch.peekLast().distance(ref.getLocus()) >= MAXIMUM_GAP_SIZE) || methyStatusListGch.size()>=MAXIMUM_DATA_POINTS){
				//send list to HMM model, get state sequence back
				if(!methyStatusListGch.isEmpty() && methyStatusListGch.size() >= MINIMUM_DATA_POINTS){
					//sendToHmm();
				//	System.err.println(ref.getLocus().getContig() + "\t" + ref.getLocus().getStart());
				}
				methyStatusListGch.clear();
				positionListGch.clear();
				//System.err.println("gapsize: " + emptyLociInList);
				emptyLociInList=0;
				return null;
			}
 		}
 		else{
 			HashSet<String> cytosinePatternConfirmedList = bvc.getSummaryAcrossRG().cytosinePatternConfirmedSet;
	 		boolean isGch = false;
	 		for(String cytosinePattern : cytosinePatternConfirmedList){
	 			if(cytosinePattern.equalsIgnoreCase("GCH")){
	 				isGch=true;
	 			}
	 			
	 		}
	 		if(isGch && bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT >= minCTdepth){
	 			//addContextToList(bvc,methyStatusListGch);
	 			double methy = (Math.random()-0.5)/100000 + (double)bvc.getSummaryAcrossRG().numC/(double)(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT);
	 			
	 			methyStatusListGch.addLast(new ObservationReal(methy));
	 			positionListGch.addLast(bvc.ref.getLocus());
	 			emptyLociInList=0;
	 		}
	 		else{
	 			emptyLociInList++;
				if(emptyLociInList >= MAXIMUM_GAP_SIZE || ((!positionListGch.isEmpty()) && positionListGch.peekLast().distance(ref.getLocus()) >= MAXIMUM_GAP_SIZE) || methyStatusListGch.size()>=MAXIMUM_DATA_POINTS){
					//send list to HMM model, get state sequence back
					if(!methyStatusListGch.isEmpty() && methyStatusListGch.size() >= MINIMUM_DATA_POINTS){
						//sendToHmm();
					//	System.err.println(ref.getLocus().getContig() + "\t" + ref.getLocus().getStart());
					}
					methyStatusListGch.clear();
					positionListGch.clear();
				//	System.err.println("gapsize: " + emptyLociInList);
					emptyLociInList=0;
					return null;
				}
	 		}
	 			
 		}
 		
		return null;
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
	public Long reduce(BisulfiteVariantCallContext value, Long sum) {
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
	
	public void onTraversalDone(Long result) {
/* Build a HMM and generate observation sequences using this HMM */
		sendToHmm();
		/*
		Hmm<ObservationReal> hmm = buildHmm();
		
		List<List<ObservationReal>> sequences;
		sequences = generateSequences(hmm);
/*
		for(List<ObservationReal> tmp : storage){
					
			System.out.println(tmp.size() + "\t" + tmp.get(0).value);
			for(ObservationReal value : tmp){
				System.out.print(value.value + ",");
			}
			System.out.println(sequences.size());
		}
		*/

		/* Baum-Welch learning */
		
		BaumWelchScaledLearner bwl = new BaumWelchScaledLearner();
		System.err.println("size:" + storage.size());
		Hmm<ObservationReal> learntHmm = buildInitHmm(storage);
		
		Hmm<ObservationReal> prevHmm = null;
		try {
			prevHmm = learntHmm.clone();
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
		while(distance >= TOLERENCE){
			i++;
			
			try {
				prevHmm = learntHmm.clone();
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			System.out.println("HMM now:\n" + prevHmm);
			
			learntHmm = bwl.iterate(learntHmm, storage);
			distance = klc.distance(learntHmm, prevHmm);
			System.out.println("Distance at iteration " + i + ": " +
					distance);
		}
		
		System.out.println("Resulting HMM:\n" + learntHmm);
		
		/* Computing the probability of a sequence */
		
		ObservationReal packetOk = new ObservationReal(0.0);
		ObservationReal packetLoss = new ObservationReal(0.5);
		
		List<ObservationReal> testSequence = 
			new ArrayList<ObservationReal>(); 
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetLoss);
		testSequence.add(packetLoss);
		testSequence.add(packetLoss);
		testSequence.add(packetLoss);
		testSequence.add(packetLoss);
		testSequence.add(packetLoss);
		testSequence.add(packetLoss);
		testSequence.add(packetLoss);
		testSequence.add(packetLoss);
		testSequence.add(packetLoss);
		testSequence.add(packetLoss);
		testSequence.add(packetLoss);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		testSequence.add(packetOk);
		
		System.out.println("Sequence probability: " +
				learntHmm.probability(testSequence));
		int[] hiddenState = learntHmm.mostLikelyStateSequence(testSequence);
		System.out.println("Sequence hidden state: ");
		for(int value : hiddenState){
			System.out.print(value);
		}
		
		System.out.println();
		
		/* Write the final result to a 'dot' (graphviz) file. */
		
		try {
			(new GenericHmmDrawerDot()).write(learntHmm, "learntHmm.dot");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
		System.out.println(storage.size());
	//	for(List<ObservationReal> tmp : storage){
	//		System.out.println(tmp.size() + "\t" + tmp.get(0).value);
	//	}
		
		/* Baum-Welch learning */
	/*	
		BaumWelchLearner bwl = new BaumWelchLearner();
		
		
		
		Hmm<ObservationReal> learntHmm = buildInitHmm();
		
		Hmm<ObservationReal> prevHmm = null;
		try {
			prevHmm = learntHmm.clone();
		} catch (CloneNotSupportedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		// This object measures the distance between two HMMs
		KullbackLeiblerDistanceCalculator klc = 
			new KullbackLeiblerDistanceCalculator();
		
		// Incrementally improve the solution
		for (int i = 0; i < 10; i++) {
			try {
				prevHmm = learntHmm.clone();
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			System.out.println("Distance at iteration " + i + ": \n" + learntHmm);
			learntHmm = bwl.iterate(learntHmm, storage);
			//System.out.println("Distance at iteration " + i + ": \n" + learntHmm);
			//		klc.distance(learntHmm, prevHmm));
		}
		
		System.out.println("Resulting HMM:\n" + learntHmm);
		
		
		//EM to estimate mixed Guassian model parameters. methylationProb
		//x = NormalDistributionImpl(double mean, double sd) 
		//x.density(methylationProb), x.cumulativeProbability(double x) 
		
		*/
		
				
				
	}
	
	//convert continous methylation level to discrete 10 different emitted states
	//Todo: change it to Guassian distribution or Beta-Binomial distribution
	private void addContextToList(BisulfiteVariantCallContext bvc, LinkedList<Integer> list){
		double methy = (double)bvc.getSummaryAcrossRG().numC/(double)(bvc.getSummaryAcrossRG().numC + bvc.getSummaryAcrossRG().numT);
		if(methy <= 0.1){
			list.addLast(0);
		}
		else if(methy > 0.1 && methy <= 0.2){
			list.addLast(1);
		}
		else if(methy > 0.2 && methy <= 0.3){
			list.addLast(2);
		}
		else if(methy > 0.3 && methy <= 0.4){
			list.addLast(3);
		}
		else if(methy > 0.4 && methy <= 0.5){
			list.addLast(4);	
		}
		else if(methy > 0.5 && methy <= 0.6){
			list.addLast(5);	
		}
		else if(methy > 0.6 && methy <= 0.7){
			list.addLast(6);
		}
		else if(methy > 0.7 && methy <= 0.8){
			list.addLast(7);
		}
		else if(methy > 0.8 && methy <= 0.9){
			list.addLast(8);
		}
		else{
			list.addLast(9);
		}
		
	}
	//to do: make whole chromosome/select CGI_TSS for the parameter training, then  input parameter as the 2nd round calculation?
	private void sendToHmm(){
		
		
		//int[][] valuesTmp = new int[1][methyStatusListGch.size()];
		//valuesTmp[numSeq] = values;
		storage.add((List<ObservationReal>) methyStatusListGch.clone());
		if(methyStatusListGch.size() > maxLen)
			maxLen = methyStatusListGch.size();
		numSeq++;
		
	}
	
/* Initial guess for the Baum-Welch algorithm */
	
	private static Hmm<ObservationReal> buildInitHmm(ArrayList<List<ObservationReal>> seqs)
	{	
		
		
		KMeansLearner<ObservationReal> kl = new KMeansLearner<ObservationReal>(STATES, new OpdfGaussianFactory(),
				seqs);
		Hmm<ObservationReal> hmm = kl.learn();
		
		
		/*
		Hmm<ObservationReal> hmm = 
			new Hmm<ObservationReal>(2,new OpdfGaussianMixtureFactory(2));
		double init1 = Math.random();
		hmm.setPi(0, init1);
		hmm.setPi(1, 1-init1);
		
		hmm.setOpdf(0, new OpdfGaussianMixture(new double[] { 0.9, 0.1 },new double[] { 1, 1 },new double[] { 0.5, 0.5 }));
		hmm.setOpdf(1,new OpdfGaussianMixture(new double[] { 0.2, 0.8 },new double[] { 1, 1 },new double[] { 0.5, 0.5 }));
		
		double init2 = Math.random();
		double init3 = Math.random();
		hmm.setAij(0, 1, init2);
		hmm.setAij(0, 0, 1-init2);
		hmm.setAij(1, 0, init3);
		hmm.setAij(1, 1, 1-init3);
		*/
		return hmm;
	}
	
	
	
	
	
	
	
	

	/* The HMM this example is based on */
	/*
	static Hmm<ObservationReal> buildHmm()
	{	
		Hmm<ObservationReal> hmm = 
			new Hmm<ObservationReal>(2,
					new OpdfGaussianMixtureFactory(2));
		
		hmm.setPi(0, 0.95);
		hmm.setPi(1, 0.05);
		
		hmm.setOpdf(0, new OpdfGaussianMixture(new double[] { 0.9, 0.1 },new double[] { 1, 1 },new double[] { 0.5, 0.5 }));
		hmm.setOpdf(1,new OpdfGaussianMixture(new double[] { 0.2, 0.8 },new double[] { 1, 1 },new double[] { 0.5, 0.5 }));
		
		hmm.setAij(0, 1, 0.05);
		hmm.setAij(0, 0, 0.95);
		hmm.setAij(1, 0, 0.10);
		hmm.setAij(1, 1, 0.90);
		
		return hmm;
	}
	*/
	
	/* Initial guess for the Baum-Welch algorithm */
	/*
	static Hmm<ObservationReal> buildInitHmm()
	{	
		Hmm<ObservationReal> hmm = 
			new Hmm<ObservationReal>(2,
					new OpdfGaussianMixtureFactory(2));
		
		hmm.setPi(0, 0.50);
		hmm.setPi(1, 0.50);
		
		hmm.setOpdf(0, new OpdfGaussianMixture(new double[] { 0.7, 0.3 },new double[] { 1, 1 },new double[] { 0.5, 0.5 }));
		hmm.setOpdf(1,new OpdfGaussianMixture(new double[] { 0.1, 0.9 },new double[] { 1, 1 },new double[] { 0.5, 0.5 }));
		
		hmm.setAij(0, 1, 0.2);
		hmm.setAij(0, 0, 0.8);
		hmm.setAij(1, 0, 0.2);
		hmm.setAij(1, 1, 0.8);
		
		return hmm;
	}
	 */
	
	/* Generate several observation sequences using a HMM */
	/*
	static <O extends Observation> List<List<O>> 
	generateSequences(Hmm<O> hmm)
	{
		MarkovGenerator<O> mg = new MarkovGenerator<O>(hmm);
		
		List<List<O>> sequences = new ArrayList<List<O>>();
		for (int i = 2; i < 200; i++)
			sequences.add(mg.observationSequence(i));

		return sequences;
	}
	*/
	/* Possible packet reception status */
	
	public enum Packet {
		OK, LOSS;
		
		public ObservationDiscrete<Packet> observation() {
			return new ObservationDiscrete<Packet>(this);
		}
	};
}
