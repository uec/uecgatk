/**
 * 
 */
package edu.usc.epigenome.uecgatk.hmm;

import java.util.Iterator;
import java.util.List;

import edu.usc.epigenome.uecgatk.distribution.OpdfBeta;
import edu.usc.epigenome.uecgatk.distribution.OpdfBetaBinomial;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Dec 11, 2012 4:20:48 PM
 * 
 */
public class ViterbiCalculator {
	/*
	 * The psy and delta values, as described in Rabiner and Juand classical
	 * papers.
	 */
	private double[][] delta; 
	private int[][] psy;
	private int[] stateSequence;
	private double lnProbability;
	
	
	/**
	 * Computes the most likely state sequence matching an observation
	 * sequence given an HMM.
	 *
	 * @param hmm A Hidden Markov Model;
	 * @param oseq An observations sequence.
	 */
	public <O extends Observation> 
	ViterbiCalculator(List<? extends O> oseq, Hmm<O> hmm)
	{
		if (oseq.isEmpty())
			throw new IllegalArgumentException("Invalid empty sequence");
		
		
		
		delta = new double[oseq.size()][hmm.nbStates()];
		psy = new int[oseq.size()][hmm.nbStates()];
		stateSequence = new int[oseq.size()];
		
		for (int i = 0; i < hmm.nbStates(); i++) {
			OpdfBeta ob = (OpdfBeta) hmm.getOpdf(i);
			int coverage = ((ObservationMethy)(oseq.get(0))).coverage;
			OpdfBetaBinomial obb = new OpdfBetaBinomial(coverage, ob.alpha(), ob.beta());
			delta[0][i] = -Math.log(hmm.getPi(i)) - 
			Math.log(obb.probability(((ObservationMethy)(oseq.get(0)))));
			psy[0][i] = 0;
		}
		
		Iterator<? extends O> oseqIterator = oseq.iterator();
		if (oseqIterator.hasNext())
			oseqIterator.next();
		
		int t = 1;
		while (oseqIterator.hasNext()) {
			O observation = oseqIterator.next();
			int coverage = ((ObservationMethy)(observation)).coverage;
			for (int i = 0; i < hmm.nbStates(); i++){
				OpdfBeta ob = (OpdfBeta) hmm.getOpdf(i);
				
				OpdfBetaBinomial obb = new OpdfBetaBinomial(coverage, ob.alpha(), ob.beta());
				double prob = obb.probability(((ObservationMethy)(observation)));
				computeStep(hmm, observation, t, i, coverage, prob);
			}
				
			
			t++;
		}
		
		lnProbability = Double.MAX_VALUE;
		for (int i = 0; i < hmm.nbStates(); i++) {
			double thisProbability = delta[oseq.size()-1][i];
			
			if (lnProbability > thisProbability) {
				lnProbability = thisProbability;
				stateSequence[oseq.size() - 1] = i;
			}
		}
		lnProbability = -lnProbability;
		
		for (int t2 = oseq.size() - 2; t2 >= 0; t2--)
			stateSequence[t2] = psy[t2+1][stateSequence[t2+1]];
	}
	
	
	/*
	 * Computes delta and psy[t][j] (t > 0) 
	 */
	private <O extends Observation> void
	computeStep(Hmm<O> hmm, O o, int t, int j, int cov, double prob) 
	{
		double minDelta = Double.MAX_VALUE;
		int min_psy = 0;
		
		for (int i = 0; i < hmm.nbStates(); i++) {
			double thisDelta = delta[t-1][i] - Math.log(hmm.getAij(i, j));
			
			if (minDelta > thisDelta) {
				minDelta = thisDelta;
				min_psy = i;
			}
		}
		
		delta[t][j] = minDelta - Math.log(prob);
		psy[t][j] = min_psy;
	}
	
	
	/**
	 * Returns the neperian logarithm of the probability of the given
	 * observation sequence on the most likely state sequence of the given
	 * HMM.
	 *
	 * @return <code>ln(P[O,S|H])</code> where <code>O</code> is the given
	 *         observation sequence, <code>H</code> the given HMM and 
	 *         <code>S</code> the most likely state sequence of this observation
	 *         sequence given this HMM.
	 */
	public double lnProbability()
	{
		return lnProbability;
	}
	
	
	/**
	 * Returns a (clone of) the array containing the computed most likely
	 * state sequence.
	 *
	 * @return The state sequence; the i-th value of the array is the index
	 *         of the i-th state of the state sequence.
	 */
	public int[] stateSequence() 
	{
		return stateSequence.clone();
	}
}
