package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import java.util.HashMap;

import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.NonRefDependSNPGenotypeLikelihoodsCalculationModel.MethylSNPModel;

public class CytosineTypeStatus {
	boolean isCpg = false;
	boolean isChh = false;
	boolean isChg = false;
	boolean isGch = false;
	boolean isGcg = false;
	boolean isHcg = false;
	double cpgMethyLevel = 0;
	double chgMethyLevel = 0;
	double chhMethyLevel = 0;
	double gchMethyLevel = 0;
	double gcgMethyLevel = 0;
	double hcgMethyLevel = 0;
	
	BisulfiteArgumentCollection BAC = null;
	
	int maxCytosineLength = 3;
	
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
		Double[] tmpDouble = new Double[3];//tmpDouble[0]: log10 likelihood in positive strand, tmpDouble[1]: log10 likelihood in negative strand, tmpDouble[2]: methylation level
		tmpDouble[0] = Double.NEGATIVE_INFINITY;
		tmpDouble[1] = Double.NEGATIVE_INFINITY;
		tmpDouble[2] = BAC.forceCpg;
		cytosineListMap.put("CG-1".toUpperCase(), tmpDouble);
		tmpDouble = new Double[3];
		tmpDouble[0] = Double.NEGATIVE_INFINITY;
		tmpDouble[1] = Double.NEGATIVE_INFINITY;
		tmpDouble[2] = BAC.forceChg;
		cytosineListMap.put("Chg-1".toUpperCase(), tmpDouble);
		tmpDouble = new Double[3];
		tmpDouble[0] = Double.NEGATIVE_INFINITY;
		tmpDouble[1] = Double.NEGATIVE_INFINITY;
		tmpDouble[2] = BAC.forceChh;
		cytosineListMap.put("Chh-1".toUpperCase(), tmpDouble);
		if(BAC.sequencingMode == MethylSNPModel.GM){
			tmpDouble = new Double[3];
			tmpDouble[0] = Double.NEGATIVE_INFINITY;
			tmpDouble[1] = Double.NEGATIVE_INFINITY;
			tmpDouble[2] = BAC.forceGch;
			cytosineListMap.put("Gch-2".toUpperCase(), tmpDouble);
			tmpDouble = new Double[3];
			tmpDouble[0] = Double.NEGATIVE_INFINITY;
			tmpDouble[1] = Double.NEGATIVE_INFINITY;
			tmpDouble[2] = BAC.forceGcg;
			cytosineListMap.put("Gcg-2".toUpperCase(), tmpDouble);
			tmpDouble = new Double[3];
			tmpDouble[0] = Double.NEGATIVE_INFINITY;
			tmpDouble[1] = Double.NEGATIVE_INFINITY;
			tmpDouble[2] = BAC.forceHcg;
			cytosineListMap.put("Hcg-2".toUpperCase(), tmpDouble);
		}
		if(!BAC.autoEstimateOtherCytosine.isEmpty()){
			String[] tmpArray = BAC.autoEstimateOtherCytosine.split(";");
			for(String tmp : tmpArray){
				tmpDouble = new Double[3];
				tmpDouble[0] = Double.NEGATIVE_INFINITY;
				tmpDouble[1] = Double.NEGATIVE_INFINITY;
				tmpDouble[2] = Double.NaN;
				cytosineListMap.put(tmp.toUpperCase(), tmpDouble);
				int tmpLength = tmp.split("-")[0].length();
				if(tmpLength > maxCytosineLength){
					maxCytosineLength = tmpLength;
				}
			}
		}
		
		if(!BAC.forceOtherCytosine.isEmpty()){
			String[] tmpArray = BAC.forceOtherCytosine.split(";");
			for(String tmp : tmpArray){
				String[] mapElement = tmp.split(":");
				tmpDouble = new Double[3];
				tmpDouble[0] = Double.NEGATIVE_INFINITY;
				tmpDouble[1] = Double.NEGATIVE_INFINITY;
				tmpDouble[2] = Double.parseDouble(mapElement[1]);
				cytosineListMap.put(mapElement[0].toUpperCase(), tmpDouble);
				int tmpLength = mapElement[0].split("-")[0].length();
				if(tmpLength > maxCytosineLength){
					maxCytosineLength = tmpLength;
				}
			}
		}
	}
	
	public CytosineTypeStatus clone(){
		CytosineTypeStatus cts = new CytosineTypeStatus(BAC);
		cts.isCpg = this.isCpg;
		cts.isChg = this.isChg;
		cts.isChh = this.isChh;
		cts.isGch = this.isGch;
		cts.isGcg = this.isGcg;
		cts.isHcg = this.isHcg;
		
		cts.cpgMethyLevel = this.cpgMethyLevel;
		cts.chgMethyLevel = this.chgMethyLevel;
		cts.chhMethyLevel = this.chhMethyLevel;
		cts.gchMethyLevel = this.gchMethyLevel;
		cts.gcgMethyLevel = this.gcgMethyLevel;
		cts.hcgMethyLevel = this.hcgMethyLevel;
		
		
		cts.BAC = this.BAC;
		
		cts.maxCytosineLength = this.maxCytosineLength;
		cts.cytosineListMap = new HashMap<String, Double[]>();
		
		for(String key : this.cytosineListMap.keySet()){
			Double[] value = this.cytosineListMap.get(key);
			cts.cytosineListMap.put(key, value);
		}
		
		return cts;
		
	}
	
}
