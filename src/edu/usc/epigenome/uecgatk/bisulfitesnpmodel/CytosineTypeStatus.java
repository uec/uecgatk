package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import java.util.HashMap;

import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.NonRefDependSNPGenotypeLikelihoodsCalculationModel.MethylSNPModel;

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

public class CytosineTypeStatus {
	boolean isC = false;
	boolean isCpg = false;
	boolean isCph = false;
	//boolean isChh = false;
	//boolean isChg = false;
	boolean isGch = false;
	boolean isHch = false;
	boolean isWch = false;
	boolean isGcg = false;
	boolean isHcg = false;
	boolean isWcg = false;
	double cytosineMethyLevel = 0;
	double cpgMethyLevel = 0;
	double cphMethyLevel = 0;
	//double chgMethyLevel = 0;
	//double chhMethyLevel = 0;
	double gchMethyLevel = 0;
	double hchMethyLevel = 0;
	double wchMethyLevel = 0;
	double gcgMethyLevel = 0;
	double hcgMethyLevel = 0;
	double wcgMethyLevel = 0;
	
	BisulfiteArgumentCollection BAC = null;
	
	int maxCytosineLength = 1;
	
	HashMap<String, Double[]> cytosineListMap = null;
	
	public CytosineTypeStatus(BisulfiteArgumentCollection bac){ 
		
		BAC = bac;
		makeCytosineMap();
	}
	
	public void makeCytosineMap(){
		makeCytosineMap(BAC);
	}
	
	public void makeCytosineMap(BisulfiteArgumentCollection BAC){
		cytosineListMap = new HashMap<String, Double[]>();
		Double[] tmpDouble = new Double[4];//tmpDouble[0]: log10 likelihood in positive strand; tmpDouble[1]: log10 likelihood in negative strand; tmpDouble[2]: methylation level; tmpDouble[3]: is it this type of cytosine: 1-true or 0-false? 
		tmpDouble[0] = Double.NEGATIVE_INFINITY;
		tmpDouble[1] = Double.NEGATIVE_INFINITY;
		tmpDouble[2] = BAC.forceCph;
		tmpDouble[3] = 0.0;
		cytosineListMap.put("C-1".toUpperCase(), tmpDouble);
		
		if(BAC.sequencingMode == MethylSNPModel.BM){
			tmpDouble = new Double[4];
			tmpDouble[0] = Double.NEGATIVE_INFINITY;
			tmpDouble[1] = Double.NEGATIVE_INFINITY;
			tmpDouble[2] = BAC.forceCpg;
			tmpDouble[3] = 0.0;
			cytosineListMap.put("CG-1".toUpperCase(), tmpDouble);
			tmpDouble = new Double[4];
			tmpDouble[0] = Double.NEGATIVE_INFINITY;
			tmpDouble[1] = Double.NEGATIVE_INFINITY;
			tmpDouble[2] = BAC.forceCph;
			tmpDouble[3] = 0.0;
			cytosineListMap.put("CH-1".toUpperCase(), tmpDouble);
			/*
			tmpDouble = new Double[4];
			tmpDouble[0] = Double.NEGATIVE_INFINITY;
			tmpDouble[1] = Double.NEGATIVE_INFINITY;
			tmpDouble[2] = BAC.forceChg;
			tmpDouble[3] = 0.0;
			cytosineListMap.put("CHG-1".toUpperCase(), tmpDouble);
			tmpDouble = new Double[4];
			tmpDouble[0] = Double.NEGATIVE_INFINITY;
			tmpDouble[1] = Double.NEGATIVE_INFINITY;
			tmpDouble[2] = BAC.forceChh;
			tmpDouble[3] = 0.0;
			cytosineListMap.put("CHH-1".toUpperCase(), tmpDouble);
			*/
		}
		
		else if(BAC.sequencingMode == MethylSNPModel.GM){
			tmpDouble = new Double[4];
			tmpDouble[0] = Double.NEGATIVE_INFINITY;
			tmpDouble[1] = Double.NEGATIVE_INFINITY;
			tmpDouble[2] = BAC.forceGch;
			tmpDouble[3] = 0.0;
			cytosineListMap.put("GCH-2".toUpperCase(), tmpDouble);
			tmpDouble = new Double[4];
			tmpDouble[0] = Double.NEGATIVE_INFINITY;
			tmpDouble[1] = Double.NEGATIVE_INFINITY;
			tmpDouble[2] = BAC.forceGch;
			tmpDouble[3] = 0.0;
			cytosineListMap.put("HCH-2".toUpperCase(), tmpDouble);
			tmpDouble = new Double[4];
			tmpDouble[0] = Double.NEGATIVE_INFINITY;
			tmpDouble[1] = Double.NEGATIVE_INFINITY;
			tmpDouble[2] = BAC.forceGch;
			tmpDouble[3] = 0.0;
			cytosineListMap.put("WCH-2".toUpperCase(), tmpDouble);
			tmpDouble = new Double[4];
			tmpDouble[0] = Double.NEGATIVE_INFINITY;
			tmpDouble[1] = Double.NEGATIVE_INFINITY;
			tmpDouble[2] = BAC.forceGcg;
			tmpDouble[3] = 0.0;
			cytosineListMap.put("GCG-2".toUpperCase(), tmpDouble);
			tmpDouble = new Double[4];
			tmpDouble[0] = Double.NEGATIVE_INFINITY;
			tmpDouble[1] = Double.NEGATIVE_INFINITY;
			tmpDouble[2] = BAC.forceHcg;
			tmpDouble[3] = 0.0;
			cytosineListMap.put("HCG-2".toUpperCase(), tmpDouble);
			tmpDouble = new Double[4];
			tmpDouble[0] = Double.NEGATIVE_INFINITY;
			tmpDouble[1] = Double.NEGATIVE_INFINITY;
			tmpDouble[2] = BAC.forceHcg;
			tmpDouble[3] = 0.0;
			cytosineListMap.put("WCG-2".toUpperCase(), tmpDouble);
		}
		if(!BAC.autoEstimateOtherCytosine.isEmpty()){
			String[] tmpArray = BAC.autoEstimateOtherCytosine.split(";");
			for(String tmp : tmpArray){
				String[] mapElement = tmp.split(":");
				tmpDouble = new Double[4];
				tmpDouble[0] = Double.NEGATIVE_INFINITY;
				tmpDouble[1] = Double.NEGATIVE_INFINITY;
				tmpDouble[2] = Double.parseDouble(mapElement[1]);
				tmpDouble[3] = 0.0;
				cytosineListMap.put(mapElement[0].toUpperCase(), tmpDouble);
				int tmpLength = mapElement[0].split("-")[0].length();
				int cytosinePos = Integer.parseInt(mapElement[0].split("-")[1]);
				if(Math.max(tmpLength-cytosinePos, cytosinePos-1) > maxCytosineLength){
					maxCytosineLength = Math.max(tmpLength-cytosinePos, cytosinePos-1);
				}
			}
		}
		
		if(!BAC.forceOtherCytosine.isEmpty()){
			String[] tmpArray = BAC.forceOtherCytosine.split(";");
			for(String tmp : tmpArray){
				String[] mapElement = tmp.split(":");
				tmpDouble = new Double[4];
				tmpDouble[0] = Double.NEGATIVE_INFINITY;
				tmpDouble[1] = Double.NEGATIVE_INFINITY;
				tmpDouble[2] = Double.parseDouble(mapElement[1]);
				tmpDouble[3] = 0.0;
				cytosineListMap.put(mapElement[0].toUpperCase(), tmpDouble);
				int tmpLength = mapElement[0].split("-")[0].length();
				int cytosinePos = Integer.parseInt(mapElement[0].split("-")[1]);
				if(Math.max(tmpLength-cytosinePos, cytosinePos-1) > maxCytosineLength){
					maxCytosineLength = Math.max(tmpLength-cytosinePos, cytosinePos-1);
				}
			}
		}
	}
	
	public CytosineTypeStatus clone(){
		CytosineTypeStatus cts = new CytosineTypeStatus(BAC.clone());
		cts.isC = this.isC;
		cts.isCpg = this.isCpg;
		cts.isCph = this.isCph;
		//cts.isChg = this.isChg;
		//cts.isChh = this.isChh;
		cts.isGch = this.isGch;
		cts.isHch = this.isHch;
		cts.isWch = this.isWch;
		cts.isGcg = this.isGcg;
		cts.isHcg = this.isHcg;
		cts.isWcg = this.isWcg;

		cts.cytosineMethyLevel = this.cytosineMethyLevel;
		cts.cpgMethyLevel = this.cpgMethyLevel;
		cts.cphMethyLevel = this.cphMethyLevel;
		//cts.chgMethyLevel = this.chgMethyLevel;
		//cts.chhMethyLevel = this.chhMethyLevel;
		cts.gchMethyLevel = this.gchMethyLevel;
		cts.hchMethyLevel = this.hchMethyLevel;
		cts.wchMethyLevel = this.wchMethyLevel;
		cts.gcgMethyLevel = this.gcgMethyLevel;
		cts.hcgMethyLevel = this.hcgMethyLevel;
		cts.wcgMethyLevel = this.wcgMethyLevel;
		
		
		cts.BAC = this.BAC.clone();
		
		cts.maxCytosineLength = this.maxCytosineLength;
		cts.cytosineListMap = new HashMap<String, Double[]>();
		
		for(String key : this.cytosineListMap.keySet()){
			Double[] value = this.cytosineListMap.get(key).clone();
			value[0] = Double.NEGATIVE_INFINITY;
			value[1] = Double.NEGATIVE_INFINITY;
			value[3] = 0.0;
			cts.cytosineListMap.put(key, value);
			
		}
		
		return cts;
		
	}
	
}
