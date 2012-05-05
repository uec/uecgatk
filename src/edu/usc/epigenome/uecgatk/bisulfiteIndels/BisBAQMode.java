/**
 * 
 */
package edu.usc.epigenome.uecgatk.bisulfiteIndels;

import java.lang.annotation.Documented;
import java.lang.annotation.ElementType;
import java.lang.annotation.Inherited;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time May 4, 2012 6:30:41 PM
 * 
 */

@Documented
@Inherited
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.TYPE)
public @interface BisBAQMode {
	public abstract edu.usc.epigenome.uecgatk.bisulfiteIndels.BisBAQ.QualityMode QualityMode() default edu.usc.epigenome.uecgatk.bisulfiteIndels.BisBAQ.QualityMode.OVERWRITE_QUALS;
    public abstract edu.usc.epigenome.uecgatk.bisulfiteIndels.BisBAQ.ApplicationTime ApplicationTime() default edu.usc.epigenome.uecgatk.bisulfiteIndels.BisBAQ.ApplicationTime.ON_INPUT;
}
