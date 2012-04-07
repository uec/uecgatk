/**
 * 
 */
package edu.usc.epigenome.uecgatk.BisSNP;

import java.util.HashMap;
import java.util.List;



/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Mar 19, 2012 9:12:14 PM
 * 
 */
public class CytosinePatternsUserDefined {

	/**
	 * 
	 */
	private List<String> cytosineContextsAcquired;
	private HashMap<String,CytosineParameters> contexts;
	private int maxCytosinePatternLen = 0;
	
	public CytosinePatternsUserDefined(List<String> cytosineContextsAcquired) {
		// TODO Auto-generated constructor stub
		this.cytosineContextsAcquired = cytosineContextsAcquired;
		contexts = new HashMap<String,CytosineParameters>();
		setContextsDefinedByUsers();
		setMaxCytosinePatternLen();
		
	}
	
	/**
	 * @return
	 */
	public HashMap<String,CytosineParameters> getContextDefined(){
		return contexts;
	}
	
	public int getMaxCytosinePatternLen(){
		return maxCytosinePatternLen;
	}
	
	public CytosinePatternsUserDefined clone(){
		return new CytosinePatternsUserDefined(cytosineContextsAcquired);
	}
	
	private void setMaxCytosinePatternLen(){

		for(String cytosinePattern : contexts.keySet()){
			int tmpLength = cytosinePattern.length();
			int cytosinePos = contexts.get(cytosinePattern).cytosinePosition;
			if(Math.max(tmpLength-cytosinePos, cytosinePos-1) > maxCytosinePatternLen){
				maxCytosinePatternLen = Math.max(tmpLength-cytosinePos, cytosinePos-1);
			}
		}
		
	}
	

	private void setContextsDefinedByUsers(){
		//System.err.println(cytosineContextsAcquired.isEmpty() + "\t" + contexts.keySet().toString());
		if(cytosineContextsAcquired.isEmpty()){
			setContextsDefinedByDefault();
		}
		
		for(String content : cytosineContextsAcquired){
			
			String[] contentsArray = content.split(",");  //-c CH,1,0.01
			CytosineParameters cp = new CytosineParameters();
			cp.cytosinePosition = Integer.parseInt(contentsArray[1]);
			if(contentsArray.length > 2){
				cp.cytosineMethylation = Double.parseDouble(contentsArray[2]);
				cp.methylationAutoestimated = false;
			}
			else{
				cp.cytosineMethylation = 0.5;
				cp.methylationAutoestimated = true;
			}
			String cytosineContext = contentsArray[0].toUpperCase();
			contexts.put(cytosineContext, cp);
			
		}

	}
	

	private void setContextsDefinedByDefault(){
		
		CytosineParameters cp1 = new CytosineParameters();
		cp1.cytosinePosition = 1;
		cp1.cytosineMethylation = 0.5;
		cp1.methylationAutoestimated = true;
		contexts.put("CG", cp1);
		CytosineParameters cp2 = new CytosineParameters();
		cp2.cytosinePosition = 1;
		cp2.cytosineMethylation = 0.5;
		cp2.methylationAutoestimated = true;
		contexts.put("CHG", cp2);
		CytosineParameters cp3 = new CytosineParameters();
		cp3.cytosinePosition = 1;
		cp3.cytosineMethylation = 0.5;
		cp3.methylationAutoestimated = true;
		contexts.put("CHH", cp3);
		CytosineParameters cp4 = new CytosineParameters();
		cp4.cytosinePosition = 1;
		cp4.cytosineMethylation = 0.5;
		cp4.methylationAutoestimated = true;
		contexts.put("CHH", cp4);
		//CytosineParameters cp3 = new CytosineParameters();
		//cp3.cytosinePosition = 1;
		//cp3.cytosineMethylation = 0.5;
		//cp3.methylationAutoestimated = true;
		//contexts.put("C", cp3);

	}
	
	//need to find a class to store all enums
	
	

}
