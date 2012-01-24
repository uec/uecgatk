package org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.CytosineTypeStatus;

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

public class BisulfiteVariantCallContext{

	public VariantContext vc = null;
    public byte refBase;
    public CytosineTypeStatus cts = null;

    // Was the site called confidently, either reference or variant?
    public boolean confidentlyCalled = false;
    public boolean emited = false;

    public BisulfiteVariantCallContext(VariantContext vc, boolean confidentlyCalledP, CytosineTypeStatus cts, boolean emited) {
        this.vc = vc;
        this.confidentlyCalled = confidentlyCalledP;
        this.cts = cts;
        this.emited = emited;
    }

    public BisulfiteVariantCallContext(VariantContext vc, byte ref, boolean confidentlyCalledP, CytosineTypeStatus cts, boolean emited) {
        this.vc = vc;
        this.refBase = ref;
        this.confidentlyCalled = confidentlyCalledP;
        this.cts = cts;
        this.emited = emited;
    }

    // blank variant context => we're a ref site
    public BisulfiteVariantCallContext(boolean confidentlyCalledP, CytosineTypeStatus cts, boolean emited) {
        this.confidentlyCalled = confidentlyCalledP;
        this.cts = cts;
        this.emited = emited;
    }

    public void setRefBase(byte ref) {
        this.refBase = ref;
    }

    public boolean isVariant() {
        if(this.vc.hasGenotypes()){
        	//System.err.println(this.vc.getGenotype(0).toString());
        	return (!this.vc.getGenotype(0).isHomRef()) && this.confidentlyCalled;
        }
        return false;
    }
    
    public boolean isHetSnp() {
    	 if(this.vc.hasGenotypes()){
         	return (this.vc.getGenotype(0).isHet()) && this.confidentlyCalled;
         }
         return false;
    }
    
}
