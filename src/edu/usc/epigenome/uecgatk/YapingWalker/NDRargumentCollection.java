package edu.usc.epigenome.uecgatk.YapingWalker;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisulfiteArgumentCollection;

public class NDRargumentCollection extends BisulfiteArgumentCollection {

	public NDRargumentCollection() {
		// TODO Auto-generated constructor stub
	}
	
	@Argument(fullName = "nucleosome_position_window", shortName = "npw", doc = "define the basic nucleosome depletion window size(bp)", required = false)
    public int nucPosWindow = 147;
	
	@Argument(fullName = "nucleosome_linker_window", shortName = "nlw", doc = "define the basic nucleosome linker window size(bp)", required = false)
    public int nucLinkerWindow = 30;
	
	@Argument(fullName = "minimum_number_gch_in_window_has_methy_value", shortName = "mgn", doc = "minimum number of gch in window has methy value", required = false)
    public int minGchNum = 5;
	
	@Argument(fullName = "minimum_CT_depth_for_gch_in_window", shortName = "mcd", doc = "minimum CT reads depth for GCH inside window", required = false)
    public int minCTDepth = 1;
	
	@Argument(fullName = "minimum_number_gch_in_linker_window_has_methy_value", shortName = "mgnlw", doc = "minimum number of gch in linker window has methy value", required = false)
    public int minGchNumLinkerWindow = 3;
	
	@Argument(fullName = "minimum_CT_depth_for_gch_in_linker_window", shortName = "mcdlw", doc = "minimum CT reads depth for GCH insidelinker  window", required = false)
    public int minCTDepthLinkerWindow = 1;
	
	@Argument(fullName = "outputFile", shortName = "outFile", doc = "bed File to which variants should be written", required = true)
    public String outFile = null;
	
	@Argument(fullName = "minimum_gch_methy_for_ndr", shortName = "ndrThreshold", doc = "minimum GCH methylation value criteria to be NDR region", required = false)
    public double ndrThreshold = 0.4;
	
	@Argument(fullName = "minimum_gch_methy_diff_for_ndr", shortName = "ndrDiffThreshold", doc = "minimum GCH methylation value differences with adjacent window to be identified as NDR region", required = false)
    public double ndrDiffThreshold = 0.4;
	
	@Argument(fullName = "enable_stat_test_for_ndr", shortName = "statTest", doc = "enable statitics test rather than hard threshold to detect NDR region", required = false)
    public boolean statTest = false;
	
	@Argument(fullName = "enable_ks_test_for_ndr", shortName = "ksTest", doc = "enable KS test rather than hard threshold to detect NDR region, otherwise, it will use Wilcoxon rank sum test", required = false)
    public boolean ksTest = false;
	
	@Argument(fullName = "sig_threshold_for_test", shortName = "sigValue", doc = "significance threshold to detect NDR region", required = false)
    public double sigValue = 0.01;
	
	@Argument(fullName = "performance_test_mode", shortName = "ptMode", doc = "enable performance test mode, which will count a bed line owns validate reads as a callable window, and output callable window also (with GCH number and CT reads depth, and average GCH methy level, CG kevele)", required = false)
    public boolean ptMode = false;
	
	@Argument(fullName = "output_callable_window_file", shortName = "ocwf", doc = "bed File name for callable window region in ptMode", required = false)
    public String ocwf = null;
	
	@Argument(fullName = "verbose_mode", shortName = "vm", doc = "enable verbose mode for debug", required = false)
    public boolean vm = false;
	
	
	public NDRargumentCollection clone() {
		NDRargumentCollection nac = new NDRargumentCollection();
		nac.nucPosWindow = nucPosWindow;
		nac.nucLinkerWindow = nucLinkerWindow;
		nac.minGchNum = minGchNum;
		nac.minCTDepth = minCTDepth;
		nac.minGchNumLinkerWindow = minGchNumLinkerWindow;
		nac.minCTDepthLinkerWindow = minCTDepthLinkerWindow;
		nac.outFile = outFile;
		nac.ndrThreshold = ndrThreshold;
		nac.ndrDiffThreshold = ndrDiffThreshold;
		nac.statTest = statTest;
		nac.ksTest = ksTest;
		nac.sigValue = sigValue;
		
		nac.ptMode = ptMode;
		nac.ocwf = ocwf;
		
		nac.vm = vm;
		return nac;
	}

}
