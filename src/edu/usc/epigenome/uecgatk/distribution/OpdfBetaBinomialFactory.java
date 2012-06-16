/**
 * 
 */
package edu.usc.epigenome.uecgatk.distribution;

import java.util.Collection;

import edu.usc.epigenome.uecgatk.NOMeSeqWalker.ObservationMethy;
import be.ac.ulg.montefiore.run.jahmm.OpdfFactory;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jun 14, 2012 7:24:51 PM
 * 
 */
public class OpdfBetaBinomialFactory implements OpdfFactory<OpdfBetaBinomial> {

	public Collection<? extends ObservationMethy> seqs;
	/**
	 * 
	 */
	public OpdfBetaBinomialFactory() {
	}
	
	//public OpdfBetaBinomialFactory(Collection<? extends ObservationMethy> co) {
	//	seqs = co;
	//}

	/* (non-Javadoc)
	 * @see be.ac.ulg.montefiore.run.jahmm.OpdfFactory#factor()
	 */
	@Override
	public OpdfBetaBinomial factor() {
		// TODO Auto-generated method stub
		return new OpdfBetaBinomial();
	}

}
