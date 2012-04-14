/**
 * 
 */
package edu.usc.epigenome.uecgatk.bisulfiteRecalibration;

import java.util.List;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.walkers.recalibration.Covariate;
import org.broadinstitute.sting.gatk.walkers.recalibration.RecalDataManager;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Apr 11, 2012 11:06:26 PM
 * 
 */
public class BisulfiteRecalDataManager extends RecalDataManager {

	/**
	 * 
	 */
	public BisulfiteRecalDataManager() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param createCollapsedTables
	 * @param numCovariates
	 */
	public BisulfiteRecalDataManager(boolean createCollapsedTables,
			int numCovariates) {
		super(createCollapsedTables, numCovariates);
		// TODO Auto-generated constructor stub
	}

	
	/**
     * Computes all requested covariates for every offset in the given read
     * by calling covariate.getValues(..).
     *
     * @param gatkRead            The read for which to compute covariate values.
     * @param requestedCovariates The list of requested covariates.
     * @return An array of covariate values where result[i][j] is the covariate
     *         value for the ith position in the read and the jth covariate in
     *         reqeustedCovariates list.
     */

    public static Comparable[][] computeCovariates(final GATKSAMRecord gatkRead, final List<Covariate> requestedCovariates, ReferenceContext ref) {
        //compute all covariates for this read
        final int numRequestedCovariates = requestedCovariates.size();
        final int readLength = gatkRead.getReadLength();

        final Comparable[][] covariateValues_offset_x_covar = new Comparable[readLength][numRequestedCovariates];
        final Comparable[] tempCovariateValuesHolder = new Comparable[readLength];

        for (int i = 0; i < numRequestedCovariates; i++) {                              // Loop through the list of requested covariates and compute the values of each covariate for all positions in this read
        	Covariate tmp = requestedCovariates.get(i);
        	if(tmp instanceof BisulfiteDinucCovariate){
        		//System.err.println(tmp.toString());
        		((BisulfiteDinucCovariate)tmp).getValues(gatkRead, tempCovariateValuesHolder,ref);
        	}
        	else{
        		tmp.getValues(gatkRead, tempCovariateValuesHolder);
        	}
        	
            for (int j = 0; j < readLength; j++)
                covariateValues_offset_x_covar[j][i] = tempCovariateValuesHolder[j];    // copy values into a 2D array that allows all covar types to be extracted at once for an offset j by doing covariateValues_offset_x_covar[j]. This avoids the need to later iterate over covar types.
        }

        return covariateValues_offset_x_covar;
    }
	
}
