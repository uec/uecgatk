package edu.usc.epigenome.uecgatk.benWalkers;

import org.broadinstitute.sting.utils.wiggle.WiggleHeader;

public class WiggleHeaderCytosines extends WiggleHeader {

	// They always make their vars private, so i have to duplicate!
    private String fName;
    // a label for the track
    private String fDescription;
    // a description of what the track is

 	public WiggleHeaderCytosines(String name, String description) {
		super(name, description);
		fName = name;
		fDescription = description;
	}

 	@Override
    public String toString() {
 		return String.format("track type=wiggle_0 name=\"%s.bare\" description=\"%s.bare\" color=204,102,0 visibility=full " +
				" graphType=points autoScale=off alwaysZero=off maxHeightPixels=64:32:10 viewLimits=-2:100", fName,fDescription);

    }

	

}
