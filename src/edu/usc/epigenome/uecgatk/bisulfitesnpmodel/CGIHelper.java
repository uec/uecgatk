package org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel;

import java.util.List;

import org.broad.tribble.Feature;

/*
 * Bis-SNP/BisSNP: It is a genotyping and methylation calling in bisulfite treated 
 * massively parallel sequencing (Bisulfite-seq and NOMe-seq) on Illumina platform
 * Copyright (C) <2011>  <Yaping Liu: lyping1986@gmail.com>

 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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


    