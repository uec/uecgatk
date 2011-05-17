package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.CytosineTypeStatus;


public class BisulfiteVariantCallContext{

	public VariantContext vc = null;
    public byte refBase;
    public CytosineTypeStatus cts = null;

    // Was the site called confidently, either reference or variant?
    public boolean confidentlyCalled = false;

    public BisulfiteVariantCallContext(VariantContext vc, boolean confidentlyCalledP, CytosineTypeStatus cts) {
        this.vc = vc;
        this.confidentlyCalled = confidentlyCalledP;
        this.cts = cts;
    }

    public BisulfiteVariantCallContext(VariantContext vc, byte ref, boolean confidentlyCalledP, CytosineTypeStatus cts) {
        this.vc = vc;
        this.refBase = ref;
        this.confidentlyCalled = confidentlyCalledP;
        this.cts = cts;
    }

    // blank variant context => we're a ref site
    public BisulfiteVariantCallContext(boolean confidentlyCalledP, CytosineTypeStatus cts) {
        this.confidentlyCalled = confidentlyCalledP;
        this.cts = cts;
    }

    public void setRefBase(byte ref) {
        this.refBase = ref;
    }

}
