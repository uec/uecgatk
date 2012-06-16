/**
 * 
 */
package edu.usc.epigenome.uecgatk.NOMeSeqWalker;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Iterator;


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
import edu.usc.epigenome.uecgatk.distribution.OpdfBetaBinomialFactory;
import edu.usc.epigenome.uecgatk.distribution.OpdfBetaFactory;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteVCFConstants;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianFactory;
import be.ac.ulg.montefiore.run.jahmm.io.FileFormatException;
import be.ac.ulg.montefiore.run.jahmm.io.HmmReader;
import be.ac.ulg.montefiore.run.jahmm.io.HmmWriter;
import be.ac.ulg.montefiore.run.jahmm.io.OpdfBetaReader;
import be.ac.ulg.montefiore.run.jahmm.io.OpdfBetaWriter;
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
	
	@Output(fullName = "states_sequence", shortName = "result", doc = "write hidden states of a sequence to a file", required = false)
	public String resultFile = null;
	
	@Argument(fullName = "min_ct_coverage", shortName = "minCT", doc = "minimum number of CT reads for count methylation level, default: 1", required = false)
	public int minCT = 1;
	
	@Argument(fullName = "training_mode", shortName = "train", doc = "enable the training mode and output HMM parameter, default: not enabled", required = false)
	public boolean train = false;
	
	@Argument(fullName = "hmm_file", shortName = "hmm", doc = "read/write HMM model from/to a file, default: read HMM parameters from a file", required = true)
	public String hmmFile = null;
	
	private static int STATES=2;
	
	private double TOLERENCE=1e-5;
	
	private int MAXIMUM_GAP_SIZE=2000000000; // maximum gap allowed for split different sequence for training & decoding
	
	private int MAXIMUM_DATA_POINTS=1000000000;
	
	private int MINIMUM_DATA_POINTS=2;
	
	private Hmm<ObservationReal> hmm = null;
	
	private bedObjectWriterImp bedWriter = null;
	
	
	
	public void initialize(){
		if(!train)
			bedWriter = new bedObjectWriterImp(new File(resultFile));
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
							Datapoint dat = new Datapoint(ref.getLocus(),new ObservationReal(methyValue));
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
		tmp.position = new ArrayList<LinkedList<GenomeLoc>>();
		tmp.value = new ArrayList<LinkedList<ObservationReal>>();
		tmp.position.add(position);
		tmp.value.add(value);
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
		if(!position.isEmpty()){
			if( position.size() >= MAXIMUM_DATA_POINTS|| position.peekLast().distance(value.position) >= MAXIMUM_GAP_SIZE || !position.peekLast().onSameContig(value.position)){
				if(position.size() >= MINIMUM_DATA_POINTS){
					sum.position.add(position);
					sum.value.add(data);
				}
				
				LinkedList<GenomeLoc> newPosition = new LinkedList<GenomeLoc>();
				LinkedList<ObservationReal> newValue = new LinkedList<ObservationReal>();
				newPosition.offerLast(value.position);
				newValue.offerLast(value.value);
				sum.position.add(newPosition);
				sum.value.add(newValue);
				return sum;
			}
		}
		position.offerLast(value.position);
		data.offerLast(value.value);
		sum.position.add(position);
		sum.value.add(data);
		
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
		for(LinkedList<ObservationReal> tmp : result.value){
			System.out.println("size: " + tmp.size());
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
			while(distance >= TOLERENCE){
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
			System.out.println("???????:\n"  );
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
			int j=0;
			for(List<ObservationReal> value : result.value){
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
				
				for(int i = 0; i < hiddenState.length; i++){
					List<Object> tmp = new LinkedList<Object>();
					tmp.add(hiddenState[i]);
					tmp.add(methyState[i]);
					bedObject bedLine = new bedObject(loci[i].getContig(), loci[i].getStart()-1, loci[i].getStop(), Strand.NONE, (List)tmp);
					bedWriter.add(bedLine);
				}
				j++;
			}

			bedWriter.close();
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
	
	
	private static Hmm<ObservationMethy> buildInitHmmByBetaBinomial(ArrayList<LinkedList<ObservationMethy>> seqs)
	{	

		KMeansLearner<ObservationMethy> kl = new KMeansLearner<ObservationMethy>(STATES, new OpdfBetaBinomialFactory(),
				seqs);
		System.out.println("KMeansLearner...");
		Hmm<ObservationMethy> hmm = kl.learn();
		return hmm;
	}
	
	public class Datum{
		public ArrayList<LinkedList<ObservationReal>> value;
		public ArrayList<LinkedList<GenomeLoc>> position;
	}
	
	public class Datapoint{
		public ObservationReal value;
		public GenomeLoc position;
		public Datapoint(GenomeLoc position, ObservationReal value){
			this.position = position;
			this.value = value;
		}
	}
	
}
