/**
 * 
 */
package edu.usc.epigenome.uecgatk.hmm;

import be.ac.ulg.montefiore.run.jahmm.ObservationReal;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Jun 15, 2012 12:01:39 PM
 * 
 */
public class ObservationMethy extends ObservationReal {

	public int coverage;
	public int distance;
	
	/**
	 * @param value
	 */
	public ObservationMethy(double value) {
		super(value);

	}
	
	/**
	 * @param value
	 */
	public void setCoverage(int coverage) {
		this.coverage = coverage;

	}

	/**
	 * @param distance to the next element
	 */
	public void setDistance(int distance) {
		this.distance = distance;

	}
}
