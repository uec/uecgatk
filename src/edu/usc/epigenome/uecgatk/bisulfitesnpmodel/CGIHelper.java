package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import java.util.List;

import org.broad.tribble.Feature;



public class CGIHelper {
	public static final String STANDARD_CGI_TRACK_NAME = "cgi";
	private CGIHelper() {}
	
	public static Feature getCGIFeature(List<Object> cgiList) {
        if (cgiList == null){
        	return null;
        }
        	
        Feature cgi = null;
        for (Object d : cgiList) {
            if (d instanceof Feature && CGIHelper.isCGI((Feature)d)) {
            	cgi = (Feature) d;
                break;
            }
        }

        return cgi;
    }
	
	public static boolean isCGI(Feature feature) {
        int size = Math.abs(feature.getEnd() - feature.getStart());
		return size > 0;
    }
	
}


    