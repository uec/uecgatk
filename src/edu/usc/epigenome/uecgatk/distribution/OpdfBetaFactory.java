/**
 * 
 */
package edu.usc.epigenome.uecgatk.distribution;

import be.ac.ulg.montefiore.run.jahmm.OpdfFactory;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jun 15, 2012 5:36:34 PM
 * 
 */
public class OpdfBetaFactory implements OpdfFactory<OpdfBeta> {

	/**
	 * 
	 */
	public OpdfBetaFactory() {
		// TODO Auto-generated constructor stub
	}

	/* (non-Javadoc)
	 * @see be.ac.ulg.montefiore.run.jahmm.OpdfFactory#factor()
	 */
	@Override
	public OpdfBeta factor() {
		return new OpdfBeta();
	}
	
	public OpdfBeta factor(int nState) {
		if(nState == 0)
			return new OpdfBeta(0.3,1.0);
		else
			return new OpdfBeta(5,1.1);
	}

}
