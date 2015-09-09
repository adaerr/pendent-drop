/*
 * Pendant drop profile integration and fitting
 *
 * This plug-in allows for interactive or automated adjustment of a
 * profile to an image of a pendant drop. The surface tension, volume
 * and surface associated with this profile can then be obtained.
 *
 * Context: The pendant drop method is commonly used to measure
 * surface tensions of liquids. It consists in analysing the shape of
 * a drop hanging typically from a capillary tube and about to detach
 * (sometimes the inverse situation of a bubble forming at the bottom
 * of a liquid is preferred, or that of a sessile drop or bubble). The
 * shape is very sensitive to the unknown interfacial tension. The
 * drop profile is described by only one non-dimensional parameter
 * (tip radius over capillary length), although in practice five
 * dimensional parameters can be adjusted within this plug-in: tip
 * position and curvature, tilt of symetry axis and capillary length.
 * The surface tension is calculated from the latter if the density
 * difference is given.
 *
 * For more information see the included PDF documentation
 * (Plugins -> Drop Analysis -> About Pendant Drop)
 */

import java.awt.Color;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.awt.geom.Point2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Area;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.HashMap;
import java.util.Iterator;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;

import net.imagej.Dataset;
import net.imagej.ImageJ;
import net.imagej.axis.Axes;
import net.imagej.axis.AxisType;
import net.imagej.display.ImageDisplay;
import net.imagej.display.OverlayService;
import net.imagej.overlay.Overlay;
import net.imagej.overlay.RectangleOverlay;
import net.imagej.table.DefaultGenericTable;
import net.imagej.table.GenericTable;

import org.scijava.ItemIO;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.command.Previewable;
import org.scijava.io.DefaultIOService;
import org.scijava.io.IOService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.service.ServiceHelper;
import org.scijava.util.Colors;
import org.scijava.util.RealRect;
import org.scijava.widget.Button;

/** An ImageJ2 plugin analyzing the shape of a pendant drop.
 */
@Plugin(type = Command.class,
        menuPath = "Plugins>Drop Analysis>Pendant Drop",
        initializer = "paramInitializer")
public class Goutte_pendante implements Command, Previewable {

    // -- Parameters --

    @Parameter
    private ImagePlus imp;

    @Parameter
    private ImageDisplay display;

    @Parameter
    private LogService log;

    @Parameter
    private OverlayService overlayService;

    @Parameter
    private StatusService statusService;

    // The 'min' attribute requires a String, but
    // Double.toString(Double.MIN_VALUE) is not a constant to the
    // compiler, so we use an explicit value close to the real
    // constant
    @Parameter(label = "Tip radius of curvature",
               persist = false, min = "1e-300")
    private double tip_radius;

    @Parameter(label = "Capillary length",
               persist = false, min = "1e-300")
    private double capillary_length;

    @Parameter(label = "Tip x coordinate",
               persist = false)
    private double tip_x;

    @Parameter(label = "Tip y coordinate",
               persist = false)
    private double tip_y;

    @Parameter(label = "Gravity angle (deg)",
               persist = false)
    private double gravity_deg;

    @Parameter(label = "Pixel size",
               initializer = "initPixelSize",
               persist = false,
               min = "1e-300")
    private double pixel_size;

    @Parameter(label = "Density contrast times g")
    private double rho_g;

    @Parameter(label = "Fit parameters checked below",
               description = "At least one parameter must be checked.",
               callback = "fitButtonCB")
    private Button fitButton;

    @Parameter(persist = false,
               visibility = org.scijava.ItemVisibility.MESSAGE,
               required = false)
    private String fittability = null;

    @Parameter(label = "Tip radius",
               callback = "fitCheckboxCB")
    private boolean fit_include_tip_radius = true;

    @Parameter(label = "Capillary length",
               callback = "fitCheckboxCB")
    private boolean fit_include_capillary_length = true;

    @Parameter(label = "Tip x coordinate",
               callback = "fitCheckboxCB")
    private boolean fit_include_tip_x = true;

    @Parameter(label = "Tip y coordinate",
               callback = "fitCheckboxCB")
    private boolean fit_include_tip_y = true;

    @Parameter(label = "Gravity angle",
               callback = "fitCheckboxCB")
    private boolean fit_include_gravity_angle = true;

    @Parameter(label = "Results", type = ItemIO.OUTPUT)
    private GenericTable results;

    // -- Other fields --

    /** A rectangle representation of dropRegion with integer precision. */
    private Rectangle bounds;
    private Area boundsArea = null;
    private Area getBoundsArea() {
        if (boundsArea == null)
            boundsArea = new Area(bounds);
        return boundsArea;
    }

    /** Arrays containing x coordinates of the drop interface. */
    private double[] leftBorder;
    private double[] rightBorder;

    // -- Constants --

    /** Enumeration of the drop shape descriptors. */
    private enum SHAPE_PARAM { TIP_RADIUS,
            CAPILLARY_LENGTH, TIP_X, TIP_Y, GRAVITY_ANGLE }

    /** Threshold defining drop interface.
     *
     * This should depend on the image type. For types other than 8bpp
     * something like (min+max)/2 could be more appropriate. TODO:
     * migrate to parameters and initialize adequately.
     */
    final double threshold = 128.f;

    /** Number of pixels around border to include in refined interface
     * detection.
     */
    final int voisinage = 10;

    /** Neighbourhood of tip for curvature estimation in dropShapeEstimator().
     */
    final int tipNeighbourhood = 5;

    /** Maximum number of iterations of Powell's minimisation
     * algorithm, in the outer loop minimizing over all parameter
     * space directions.
     */
    final int maxIterDir = 100;

    /** Maximum number of iterations inside the 1D minimisation
     * function for each direction.
     */
    final int maxIterMin = 100;

    // -- Command methods --

    @Override
    public void run() {
        HashMap<SHAPE_PARAM,Boolean> fitMe = tagParamsToFit();
        if ( ! fitMe.containsValue(Boolean.TRUE) ) {
            log.error("At least one parameter must be selected !");
            return;
        }

        Contour drop = new Contour(tip_radius, capillary_length,
                                   tip_x, tip_y, gravity_deg);

        ImageStack stack = imp.getStack();

        // create results table
        final int colCount = 10;
        results = new DefaultGenericTable(colCount, stack.getSize());
        results.setColumnHeader(0, "slice");
        results.setColumnHeader(1, "tip radius");
        results.setColumnHeader(2, "tip x");
        results.setColumnHeader(3, "tip y");
        results.setColumnHeader(4, "gravity angle");
        results.setColumnHeader(5, "capillary length");
        results.setColumnHeader(6, "surface tension");
        results.setColumnHeader(7, "volume");
        results.setColumnHeader(8, "surface");
        results.setColumnHeader(9, "fit distance");

        // iterate over slices of the stack
        for (int n=0; n<stack.getSize(); n++) {
            statusService.showStatus("Pendant drop processing image " + (n+1));

            if (! findDropBorders(stack.getProcessor(n+1))) {
                log.error("Could not detect drop on image " + (n+1) + ", skipping.");
                continue;
            }

            final ContourProperties dropFit =
                fitContourToImage(drop, fitMe, false);

            // take this contour as starting point for fitting the next image
            drop = dropFit.getContour();

            // copy parameters into results table
            results.set(0, n, n+1);
            results.set(1, n, drop.getTipRadius());
            results.set(2, n, drop.getTipX());
            results.set(3, n, drop.getTipY());
            results.set(4, n, drop.getGravityAngle());
            final double capLength = drop.getCapillaryLength();
            results.set(5, n, capLength);
            final double surface_tension = capLength * capLength * rho_g;
            results.set(6, n, surface_tension);
            results.set(7, n, dropFit.getVolume());
            results.set(8, n, dropFit.getSurface());
            results.set(9, n, dropFit.getFitDistance());

            imp.setSlice(n+1);
            updateOverlay(dropFit.getDimShape());
            imp.updateAndDraw();

            log.info("Properties of image "+(n+1)+":");
            showProperties(dropFit);
        }

        statusService.showStatus("Pendant drop done");
    }

    // -- Previewable methods --

    @Override
    public void preview() {
        Contour c = new Contour(tip_radius, capillary_length, tip_x,
                                tip_y, gravity_deg);
        log.info("previewing contour "+c);
        ContourProperties cp = new ContourProperties(c);
        updateOverlay(cp.getDimShape());
        showProperties(cp);
    }

    @Override
    public void cancel() {
        log.info("cancelled");
    }

    // -- Initializer methods --

    /** Initializes some parameters by roughly analyzing the image.
     * The corresponding parameters are: {@link #tip_radius},
     * {@link #capillary_length}, {@link #tip_x}, {@link #tip_y},
     * {@link #gravity_deg}
     */
    protected void paramInitializer() {
        // get selection bounds as rectangle
        RealRect r = overlayService.getSelectionBounds(display);
        bounds = new Rectangle((int) Math.round(r.x),
                               (int) Math.round(r.y),
                               (int) Math.round(r.width),
                               (int) Math.round(r.height));
        // work around getSelectionBounds providing old image
        // dimensions after a crop
        if (bounds.x + bounds.width > imp.getWidth() ||
            bounds.y + bounds.height > imp.getHeight()) {
            bounds.x = 0;
            bounds.y = 0;
            bounds.width = imp.getWidth();
            bounds.height = imp.getHeight();
        }
        // work around getSelectionBounds not seeing an imp's roi
        if (imp.getRoi() != null) {
            bounds = imp.getRoi().getBounds();
        }
        log.info("drop region: +" + bounds.x + " +" +  bounds.y
                 + ", " +  bounds.width + " x " + bounds.height);

        // pixel size and density contrast parameters
        ij.measure.Calibration cal = imp.getCalibration();
        if (cal.scaled()) {
            pixel_size = cal.pixelWidth;
            rho_g = 9.81;
        } else {
            pixel_size = 1;
            rho_g = 1;
        }

        // parameters caracterizing the drop shape
        dropShapeEstimator(imp.getProcessor());
    }

    /** Analyze given image to roughly estimate drop shape descriptors. */
    protected void dropShapeEstimator(ImageProcessor ip) {

        if (! findDropBorders(ip)) {
            log.error("Could not detect drop, unrecoverable plugin failure.");
            throw new RuntimeException("Could not detect drop.");
        }
        final int tip = leftBorder.length - 1;

        // tip curvature: linear fit to half width squared near tip
        final double[] radiusSquare = new double[tipNeighbourhood];
        final double[] yCoord = new double[tipNeighbourhood];
        for (int i=0; i<tipNeighbourhood; i++) {
            double halfWidth = 0.5*( rightBorder[tip-i] - leftBorder[tip-i] );
            radiusSquare[i] = halfWidth * halfWidth;
            yCoord[i] = i;
        }
        Polynome tipFit = linearFit(yCoord, radiusSquare);
        tip_radius = 0.5 * Math.abs(tipFit.getCoeff(1)) * pixel_size;

        // in case estimation fails badly assign some reasonable value
        if (Double.isNaN(tip_radius) || tip_radius <= 0 ||
            tip_radius > Math.min(bounds.width,bounds.height) * pixel_size)
            tip_radius = tip * pixel_size / 4;

        // direction of gravity: tilt of center line
        //
        // link middle of top line to estimated center of drop
        // (one radius above tip)
        double dX = 0.5*( rightBorder[0] + leftBorder[0]
                                - rightBorder[tip] - leftBorder[tip] );
        // in case estimation fails badly assign some reasonable value
        if (Double.isNaN(dX)) dX = 0;
        final double dY = tip - tip_radius / pixel_size;
        final double tiltAngleRad = Math.atan2(-dX,dY);
        gravity_deg = tiltAngleRad * 180/Math.PI;

        // drop tip
        tip_x = ( bounds.x +
                  0.5*( rightBorder[tip] + leftBorder[tip] )) * pixel_size;
        tip_y = ( bounds.y + tip ) * pixel_size;
        final double tip_correction = tipFit.getCoeff(0) / tipFit.getCoeff(1);
        if (!Double.isNaN(tip_correction) && Math.abs(tip_correction) < tipNeighbourhood)
            tip_y += tip_correction * pixel_size;
        // the former expression assumes the drop is vertical, but we can
        // correct for tilt:
        tip_x -= tip_radius * Math.sin(tiltAngleRad);
        tip_y -= tip_radius * (1 - Math.cos(tiltAngleRad));

        // in case estimation fails badly assign some reasonable value
        if (tip_x/pixel_size < bounds.x || tip_x/pixel_size >= bounds.x + bounds.width)
            tip_x = ( bounds.x + bounds.width/2 ) * pixel_size;
        if (tip_y < bounds.y || tip_y >= bounds.y + bounds.height)
            tip_y = ( bounds.y + tip ) * pixel_size;

        // capillary length: from curvature difference between tip
        // (=1/tip_radius) and the point of maximum horizontal drop
        // diameter (to be calculated here)

        final int bellyNeighbourhood = tip/8;// fit on up to twice
        // these many points

        // find maximum horizontal extent of drop
        int yBelly = tip;
        for (int i=1; i < tip; i++) {
            double hDiam = rightBorder[tip-i] - leftBorder[tip-i];
            if (Double.isNaN(hDiam)) continue;
            if (hDiam > rightBorder[yBelly] - leftBorder[yBelly])
                yBelly = tip-i;
        }
        // define neighbourhood for fit
        if (yBelly < bellyNeighbourhood)// too close to top;
            yBelly = bellyNeighbourhood;
        int yBelow = yBelly - bellyNeighbourhood;
        int yAbove = Math.min(tip, yBelly + bellyNeighbourhood);
        double[] xB = new double[yAbove - yBelow + 1];
        for (int i = yBelow; i <= yAbove; i++)
            xB[i-yBelow] = 0.5 * (rightBorder[i] - leftBorder[i]);
        double[] yB = new double[yAbove - yBelow + 1];
        for (int i = yBelow; i <= yAbove; i++)
            yB[i-yBelow] = i - yBelly;

        //log.info("tip = "+tip+", yBelly = "+yBelly);

        // fit parabola to neighbourhood of yBelly
        Polynome bellyFit = quadraticFit(yB, xB);

        // and calculate curvatures at yBelly
        double c1 = 1 / bellyFit.getValueAt(0); // horizontal curvature
        double c2 = - 2*bellyFit.getCoeff(2) /
            Math.cos(Math.atan(bellyFit.getCoeff(1)));
        double c0 = pixel_size / tip_radius;// tip curvature (pix^-1)

        //log.info("c0 = "+c0+", c1 = "+c1+", c2 = "+c2);

        // pressure difference between yBelly and tip gives capillary
        // length estimate
        capillary_length = Math.sqrt( (tip-yBelly) / (2*c0-c1-c2) ) *pixel_size;

        // in case estimation fails badly assign some reasonable value
        if (Double.isNaN(capillary_length) || capillary_length == 0)
            capillary_length = 2*tip_radius;
    }

    /** Detect drop borders and store positions in the
     * left/rightBorder arrays.
     */
    private boolean findDropBorders(ImageProcessor ip) {
        leftBorder = null;
        rightBorder = null;

        for (int y = bounds.height - 1; y >= 0; y--) {

            // find border positions with integer precision

            // left border first
            int xl = 0;
            while (xl < bounds.width &&
                   ip.getPixelValue(bounds.x + xl, bounds.y + y) > threshold)
                xl ++;

            if (xl >= bounds.width) {// drop not detected in this scanline
                if (leftBorder != null) {
                    leftBorder[y]  = Double.NaN;
                    rightBorder[y] = Double.NaN;
                }
                continue;
            } else if (leftBorder == null) {
                // allocate array on drop tip detection
                leftBorder = new double[y+1];
                rightBorder = new double[y+1];
            }

            // right border next
            int xr = bounds.width - 1;
            while (xr > xl &&
                   ip.getPixelValue(bounds.x + xr, bounds.y + y) > threshold)
                xr --;
            xr ++; // so xl and xr point just to the right of the interface

            // don't go further if not enough pixels for subpixel-fitting
            if (xr - xl <= voisinage ||
                xl - voisinage < 0 || xr + voisinage > bounds.width) {
                leftBorder[y]  = xl - 0.5;
                rightBorder[y] = xr - 0.5;
                continue;
            }

            // now determine drop borders with sub-pixel precision
            leftBorder[y]  = fitStep(ip, xl, y, voisinage, false);
            rightBorder[y] = fitStep(ip, xr, y, voisinage, true);
        } // end for y

        if (leftBorder == null)
            return false;
        else
            return true;
    }

    /** Calculate sub-pixel position of a (rising/falling) step having
     * the same integral as the image intensity profile around the
     * position (x,y). Used to refine drop interface detection.
     */
    private double fitStep(ImageProcessor ip, int x, int y, int n,
                           boolean rising) {
        double acc = 0;
        double minValue = Double.POSITIVE_INFINITY;
        double maxValue = Double.NEGATIVE_INFINITY;
        for (int dx = -n; dx < n; dx++) {
            double v = ip.getPixelValue(bounds.x + x + dx, bounds.y + y);
            acc += v;
            if (v > maxValue) maxValue = v;
            if (v < minValue) minValue = v;
        }
        acc -= 2*n*minValue;
        acc /= maxValue - minValue;
        if (rising)
            return x - 0.5 + n - acc;
        else
            return x - 0.5 - n + acc;
    }

    // -- Callbacks --

    /** Called by the button that triggers a fit of the drop contour. */
    private void fitButtonCB() {
        HashMap<SHAPE_PARAM,Boolean> fitMe = tagParamsToFit();
        if ( ! fitMe.containsValue(Boolean.TRUE) ) {
            log.error("At least one parameter must be selected !");
            return;
        }

        Contour drop = new Contour(tip_radius, capillary_length,
                                   tip_x, tip_y, gravity_deg);

        final ContourProperties dropFit =
            fitContourToImage(drop, fitMe, true);

        // It's a programming error to modify parameter values from a callback.
        // Have to do this properly one day (DynamicCommand?), but in
        // the meantime all we care is to store this (hopefully
        // improved) result as a starting point for coming
        // calculations.
        drop = dropFit.getContour();
        tip_radius = drop.getTipRadius();
        capillary_length = drop.getCapillaryLength();
        tip_x = drop.getTipX();
        tip_y = drop.getTipY();
        gravity_deg = drop.getGravityAngle();

        updateOverlay(dropFit.getDimShape());
        imp.updateAndDraw();
        showProperties(dropFit);
    }

    /** Check if at least one fit parameter is checked.
     */
    private void fitCheckboxCB() {
        if (!( fit_include_tip_radius || fit_include_capillary_length ||
               fit_include_tip_x || fit_include_tip_y ||
               fit_include_gravity_angle )) {
            if (fittability == null)
                fittability = "At least one parameter must be checked !";
        } else {
            if (fittability != null) fittability = null;
        }
    }

    // -- Processing --

    /** Build a map from the input parameters, indicating for each
     * shape parameter whether it should be adjusted during fitting
     * (true) or left unchanged (false). Such a map is required as
     * input parameter to the fitContourToImage() method.
     *
     * @see SHAPE_PARAM
     * @see fit_include_tip_radius
     * @see fit_include_capillary_length
     * @see fit_include_tip_x
     * @see fit_include_tip_y
     * @see fit_include_gravity_angle
     * @see fitContourToImage()
     */
    private HashMap<SHAPE_PARAM,Boolean> tagParamsToFit() {
        final HashMap<SHAPE_PARAM,Boolean> fitMe =
            new HashMap<SHAPE_PARAM,Boolean>();
        fitMe.put(SHAPE_PARAM.TIP_RADIUS, fit_include_tip_radius);
        fitMe.put(SHAPE_PARAM.CAPILLARY_LENGTH, fit_include_capillary_length);
        fitMe.put(SHAPE_PARAM.TIP_X, fit_include_tip_x);
        fitMe.put(SHAPE_PARAM.TIP_Y, fit_include_tip_y);
        fitMe.put(SHAPE_PARAM.GRAVITY_ANGLE, fit_include_gravity_angle);
        return fitMe;
    }

    /** Hydrostatic pressure equilibrium equations for an axisymetric drop.
     * variables (r, psi, kappa, z) represent
     *
     * r ... drop radius in cylindrical coordinates
     * psi ... interface tilt to horizontal
     * kappa ... radial curvature
     * z ... height
     */
    private void deriv(double state[], double deriv[]) {
        final double c = Math.cos(state[1]);
        final double s = Math.sin(state[1]);
        deriv[0] = c;
        deriv[1] = state[2];
        deriv[2] = - s + c*(s/state[0] - state[2])/state[0];
        deriv[3] = s;
    }

    /** Calculate non-dimensional drop profile through integration of
     * the hydrostatic pressure equilibrium. Assumes axisymetry.
     * Expects tip radius and maximum z value to be given normalised
     * by capillary length. */
    ArrayList<Point2D.Double> calculateProfile(double tipRadius, double zMax) {
        //System.err.println("tipRadius = "+tipRadius+", zMax = "+zMax);

        final double ds = 0.001;// integration step
        final ArrayList<Point2D.Double> profile = new ArrayList<Point2D.Double>((int)(zMax/ds));

        // Near the drop tip the integral equations become singular,
        // but the solution is simply a nearly perfectly spherical
        // profile (at least at distances small wrt capillary length
        // and tip radius): we use an approximate solution to get away
        // from the singularity at r=0.
        //
        // The following constant defines the fraction of capillary
        // length or tip radius (whichever is smaller) up to which we
        // use the approximate solution
        final double approxLimit = 1;//0.01;

        double s;
        double sLimit;// curvilinear coordinate where we switch from
        // the approximate solution to integration
        if (tipRadius < 1)
            sLimit = approxLimit*tipRadius;// tipRadius is limiting
        else
            sLimit = approxLimit;// the capillary length is limiting
        profile.add(new Point2D.Double(0, 0));
        for (s=0; s < sLimit; s += ds) {
            double r, z;
            r = tipRadius*Math.sin(s/tipRadius) +
                Math.pow(s,5)/(40*tipRadius*tipRadius);
            z = 0.5f*s*s*(1-s*s*(0.75f+1/(tipRadius*tipRadius))/12)/tipRadius;
            profile.add(new Point2D.Double((double)r,(double)z));
        }

        // now continue through integration with a Runge-Kutta scheme
        double[] state, tmp, k1, k2, k3, k4;
        int nVar = 4;// r, psi, kappa, z
        state = new double[nVar];
        tmp = new double[nVar];
        k1 = new double[nVar];
        k2 = new double[nVar];
        k3 = new double[nVar];
        k4 = new double[nVar];
        // intial conditions using approximate solution
        state[0] = tipRadius*Math.sin(s/tipRadius) +
            Math.pow(s,5)/(40*tipRadius*tipRadius);// r
        state[1] = s*(1 - 0.125f*s*s)/tipRadius;// psi
        state[2] = (1 - 0.375f*s*s)/tipRadius;// kappa
        state[3] = 0.5f*s*s*(1 - s*s*(0.75f + 1/(tipRadius*tipRadius))/12) /
            tipRadius;// z
        while (state[3] < zMax) {// integrate until tall enough
            deriv(state,k1);
            for (int i=0; i<nVar; i++) tmp[i] = state[i] + 0.5f*ds*k1[i];
            deriv(tmp,k2);
            for (int i=0; i<nVar; i++) tmp[i] = state[i] + 0.5f*ds*k2[i];
            deriv(tmp,k3);
            for (int i=0; i<nVar; i++) tmp[i] = state[i] + ds*k3[i];
            deriv(tmp,k4);
            final boolean thinning = 2*state[1] > Math.PI;
            for (int i=0; i<nVar; i++)
                state[i] += ds*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
            final boolean hasNeck = thinning && (2*state[1] < Math.PI);
            final double dz = ds*(k1[3]+2*k2[3]+2*k3[3]+k4[3])/6;
            // abort if profile becomes non-physical:
            if (dz < 0 || //  bends back downwards
                state[0] <= 0 || // traverses axis 
                hasNeck) {
                break;
            }
            profile.add(new Point2D.Double((double)state[0],(double)state[3]));
        }

        return profile;
    }

    /** Construct closed drop contour from profile by concatenating
     * the profile with its reflection about the y-axis
     */
    private Path2D makeClosedPath(final ArrayList<Point2D.Double> profile) {
        final Path2D.Double drop = new Path2D.Double();
        Point2D.Double p;
        int N = profile.size();
        p = profile.get(0);
        drop.moveTo(p.x, p.y);
        for (int i = 1; i < N; i++) {
            p = profile.get(i);
            drop.lineTo(p.x, p.y);
        }
        for (int i = N-1; i > 0; i--) {
            p = profile.get(i);
            drop.lineTo(-p.x, p.y);
        }
        drop.closePath();

        return drop;
    }

    /** Transform the non-dimensional drop path into screen space.
     *
     *  @param drop non-dimensional drop path
     *  @param scale capillary length in pixels
     *  @param tipX tip x coordinate in pixels
     *  @param tipY tip y coordinate in pixels
     *  @param angle_deg deviation of gravity from vertical in degrees
     *  @param clipArea intersect resulting shape with this area
     */
    Shape contourToScreen(Path2D drop, double scale, double tipX, double tipY,
                          double angle_deg, Area clipArea) {

        // read this bottom-up to get the order of transformations right...
        AffineTransform toScreenT;
        toScreenT = AffineTransform.getTranslateInstance(tipX, tipY);
        toScreenT.scale(scale, scale);
        toScreenT.rotate(angle_deg * Math.PI/180);
        toScreenT.quadrantRotate(2);
        //return drop.createTransformedShape(toScreenT);

        // Generate a low resolution version of the drop path. If the
        // path is polygonal with roughly pixel precision we will not
        // see the difference in the fit (except on the computation
        // time...)
        final double flatness = 0.5; // Probably no effect unless the
                                     // flattening path iterator also
                                     // flattens polygons.
        final double minDist2 = 1.0; // Square of minimal distance
                                     // between points in the low res
                                     // version (in pixels).
        Path2D lowresDrop = new Path2D.Double();
        double[] coords = new double[6];
        double ox, oy;
        PathIterator iter = drop.getPathIterator(toScreenT, flatness);
        iter.currentSegment(coords);
        ox = coords[0];
        oy = coords[1];
        lowresDrop.moveTo(ox, oy);
        iter.next();
        while (!iter.isDone()) {
            int type = iter.currentSegment(coords);
            if (type == PathIterator.SEG_CLOSE) {
                lowresDrop.closePath();
                break;
            }
            double x = coords[0];
            double y = coords[1];
            //System.err.println("type = "+type+", x = "+x+", y = "+y+
            //       ", d2 = "+((x-ox)*(x-ox)+(y-oy)*(y-oy)));
            if ((x-ox)*(x-ox)+(y-oy)*(y-oy) > minDist2) {
                lowresDrop.lineTo(x,y);
                ox = x;
                oy = y;
            }
            iter.next();
        }
        lowresDrop.closePath();

        Area visibleDrop = new Area(lowresDrop);
        visibleDrop.intersect(clipArea);

        return visibleDrop;
    }

    /** Calculate drop contour volume (non-dimensional). */
    private double calcVolume(final ArrayList<Point2D.Double> profile,
                              double yMax) {
        Point2D.Double p, q;
        double volume = 0;
        Iterator<Point2D.Double> iter = profile.iterator();
        p = iter.next();
        while (iter.hasNext()) {
            q = iter.next();
            if (q.y > yMax) {log.info("q.y=" +q.y+ ", yMax = " +yMax);break;}
            volume += (Math.PI/3) * (p.x*p.x + p.x*q.x + q.x*q.x) * (q.y - p.y);
            p = q;
        }
        return volume;
    }

    /** Calculate drop contour surface (non-dimensional). */
    private double calcSurface(final ArrayList<Point2D.Double> profile,
                               double yMax) {
        Point2D.Double p, q;
        double surface = 0;
        Iterator<Point2D.Double> iter = profile.iterator();
        p = iter.next();
        while (iter.hasNext()) {
            q = iter.next();
            if (q.y > yMax) break;
            surface += Math.PI * (p.x + q.x) *
                Math.sqrt((q.x - p.x) * (q.x - p.x) +
                          (q.y - p.y) * (q.y - p.y));
            p = q;
        }
        return surface;
    }

    double calcFitDistance(Shape drop) {
        double[] leftIntersect = new double[bounds.height];
        double[] rightIntersect = new double[bounds.height];
        for (int y = 0; y < bounds.height; y++) {
            leftIntersect[y] = Double.NaN;
            rightIntersect[y] = Double.NaN;
        }
        // follow the contour and note intersection point for every line
        PathIterator pi = drop.getPathIterator(null);
        if (pi.isDone()) return Double.POSITIVE_INFINITY;
        double[] coord = new double[6];
        double xo, yo, xn, yn;
        int segType, yof, ynf;
        segType = pi.currentSegment(coord);
        xo = coord[0] - bounds.x;
        yo = coord[1] - bounds.y;
        yof = (int)Math.floor(yo);
        leftIntersect[(int)yo] = xo;

        pi.next();

        while (!pi.isDone()) {
            segType = pi.currentSegment(coord);
            if (segType != PathIterator.SEG_LINETO) {
                pi.next();
                continue;
            }
            xn = coord[0] - bounds.x;
            yn = coord[1] - bounds.y;
            // find out if this segment crossed a horizontal
            //line of integer coordinate
            ynf = (int)Math.floor(yn);
            if (ynf != yof) {// integer bound crossed !
                int ymin = Math.max(0, Math.min(ynf, yof) + 1);
                int ymax = Math.min(bounds.height - 1, Math.max(ynf, yof));
                for (int yi = ymin; yi <= ymax; yi++) {
                    double frac = (yi - yo) / (yn - yo);
                    double xi = xo + frac * (xn - xo);
                    if (ynf > yof) // left side
                        leftIntersect[yi] = xi;
                    else
                        rightIntersect[yi] = xi;
                }
            }
            //log.info("lineto (" + xn + ", " + yn + ", " + ynf + ")");
            xo = xn;
            yo = yn;
            yof = ynf;
            pi.next();
        }
        rightIntersect[(int)yo] = xo;

        // add values where there is no contour
        double dummyIntersect = bounds.width/2;
        for (int y=0; y < bounds.height; y++) {
            if (Double.isNaN(leftIntersect[y]) || Double.isNaN(rightIntersect[y])) {
                if (Double.isNaN(leftIntersect[y]))
                    leftIntersect[y] = dummyIntersect;
                if (Double.isNaN(rightIntersect[y]))
                    rightIntersect[y] = dummyIntersect;
            } else {
                dummyIntersect = 0.5*(leftIntersect[y] + rightIntersect[y]);
            }
            //log.info("y = " + y + ": left = " + leftIntersect[y] + ", right = " + rightIntersect[y]);
        }

        // calculate penalty by summing over the distances between the measured
        // and the theoretical border position.
        double ChiSq = 0;
        double xml = 0;
        double xtl = 0;
        double xmr = 0;
        double xtr = 0;
        for (int y = 0; y < leftBorder.length; y++){
          if (!Double.isNaN (leftBorder[y])){
              xml = leftBorder[y];
            if (!Double.isNaN (leftIntersect[y])){
                xtl = leftIntersect[y];
            } else {
                xtl = bounds.width/2;
            }
            ChiSq += (xml - xtl)*(xml - xtl); //alternative: abs(xml - xtl)
           }

        if (!Double.isNaN(rightBorder[y])){
              xmr = rightBorder[y];
              if (!Double.isNaN(rightIntersect[y])){
                  xtr = rightIntersect[y];
              } else {
                  xtr = bounds.width/2;
              }
            ChiSq += (xmr - xtr)*(xmr - xtr); //alternative: abs(xmr - xtr)
            }
        }
        return ChiSq;
    }

    /** Set given curve as overlay on image. */
    private void updateOverlay(Shape c) {
        //log.info("updating overlay");
        if (c == null) return;
        // translation by +(0.5,0.5)
        AffineTransform t = new AffineTransform(1, 0, 0, 1, 0.5, 0.5);
        //imp.setOverlay(t.createTransformedShape(c), Color.red, null);
        Path2D.Float border = new Path2D.Float();
        border.moveTo(bounds.x + leftBorder[0], bounds.y + 0);
        for (int y = 0; y < leftBorder.length; y++) {
        border.lineTo(bounds.x + leftBorder[y], bounds.y + y);
        }
        for (int y = rightBorder.length-1; y>=0; y--) {
        border.lineTo(bounds.x + rightBorder[y], bounds.y + y);
        }
        border.closePath();

        ij.gui.Roi r = new ij.gui.ShapeRoi(t.createTransformedShape(new Area(border)));
        r.setStrokeColor(Color.blue);
        ij.gui.Overlay o = new ij.gui.Overlay(r);

        r = new ij.gui.ShapeRoi(t.createTransformedShape(c));
        r.setStrokeColor(Color.red);
        o.add(r);

        imp.setOverlay(o);
    }

    private void showProperties(ContourProperties drop) {
        final double capillary_length = drop.getContour().getCapillaryLength();
        final double surface_tension = capillary_length*capillary_length* rho_g;
        log.info("surface tension = " + surface_tension +
                 "\ndrop volume = " + drop.getVolume() +
                 "\ndrop surface = " + drop.getSurface() +
                 "\nfitDistance = " + drop.getFitDistance());
    }

    /** Find a drop shape that minimises the distance to the given
     * image, by varying a subset of the shape parameters. This is
     * done by iteratively 1D-minimising along selected directions,
     * where the latter are updated according to a method described by
     * Powell.
     *
     * @param cd a drop shape used as the starting point
     * @param ip the image of the measured drop
     * @param fitMe a map indicating which shape parameters to fit
     * @param drawOverlay fit silently (false) or update profile continuously (true)
     * @see tagParamsToFit()
     */
    public ContourProperties fitContourToImage(Contour cd,
                                               HashMap<SHAPE_PARAM,Boolean> fitMe,
                                               boolean drawOverlay) {

        final double precision1 = 1e-8;
        final double precision2 = 1e-15;

        Contour xCur, xOld, xE;
        ContourProperties pCur, pOld, pE, pMin;
        double fCur, fOld, fE;

        xCur = cd;
        pCur = new ContourProperties(xCur);
        fCur = pCur.getFitDistance();

        int Ndirs = 0;
        for (SHAPE_PARAM p: SHAPE_PARAM.values())
            if (fitMe.get(p)) Ndirs++;
        // create list of directions
        ArrayList<DContour> directions = new ArrayList<DContour>(Ndirs);
        for (SHAPE_PARAM p: fitMe.keySet()) {
            if (fitMe.get(p)) {
                DContour direction = new DContour(0, 0, 0, 0, 0);
                if (p == SHAPE_PARAM.GRAVITY_ANGLE &&
                    fitMe.get(SHAPE_PARAM.TIP_X)) {
                    direction.put(p, 1/cd.get(SHAPE_PARAM.TIP_RADIUS));
                    direction.put(SHAPE_PARAM.TIP_X, 1.0);
                } else {
                    direction.put(p, 1.0);
                }
                directions.add(direction);
            }
        }
        int step;
        for (step=0; step < maxIterDir; step++) {
            //log.info("----- new round -----");
            fOld = fCur;
            xOld = xCur;
            double deltaFmax = 0;
            DContour maxDir = null;// direction of steepest descent
            int numDownSteps = 0;// number of descending steps this round
            // loop over directions and find direction which yields
            // the largest function decrease
            for (DContour dir: directions) {
                pMin = minimise1D(xCur, fCur, dir);
                double deltaF = fCur - pMin.getFitDistance();
                //log.info("dir "+dir+": fCur = "+fCur+": deltaF = "+deltaF+", deltaFmax = "+deltaFmax);
                if (deltaF > deltaFmax) {
                    maxDir = dir;
                    deltaFmax = deltaF;
                }
                if (deltaF > 0) {
                    pCur = pMin;
                    xCur = pMin.getContour();
                    fCur = pMin.getFitDistance();
                    numDownSteps ++;
                }
                if (drawOverlay) updateOverlay(pCur.getDimShape());
            }

            // test for convergence
            //log.info("change = "+Math.abs(fCur - fOld)+", threshold = "+(precision1*Math.abs(fCur) + precision2));
            if (Math.abs(fCur - fOld) < precision1*Math.abs(fCur) + precision2)
                break;

            // construct extrapolated point and direction from xCur and xOld
            if (numDownSteps == 0) continue;
            DContour newDirection = new DContour(xOld, xCur);
            //log.info("testing direction: "+newDirection);
            xE = xCur.plus(newDirection);
            pE = new ContourProperties(xE);
            fE = pE.getFitDistance();
            //log.info("fOld = "+fOld+", fCur = "+fCur+", fE = "+fE);

            // now there are 3 points: (xOld, fOld), (xCur, fCur)
            // and (xE, fE); check where minimum is on this line
            if (fE < fCur) {
                double test = 2*(fOld - 2*fCur + fE) *
                    (fOld - fCur - deltaFmax) * (fOld - fCur - deltaFmax) -
                    (fOld - fE) * (fOld - fE) * deltaFmax;
                //log.info("test = "+test);
                if (test < 0.0) {// use new direction
                    //log.info("new direction accepted, removing direction " +
                    //         maxDir);
                    xCur = xE;
                    pCur = pE;
                    fCur = fE;
                    directions.remove(maxDir);
                    directions.add(0, newDirection);// prepend to list
                }
            }
        }
        if (step >= maxIterDir)
            log.error("minimisation failed");
        //else
        //    log.info("minimisation succeeded");
        return pCur;
    }

    /** Calculate where the minimum of the parabola going through the
     * three given points is. No checks are performed to make sure the
     * input points are not singular (e.g. x1==x2). */
    private double minParabola(double x1, double x2, double x3,
                       double y1, double y2, double y3) {
        double x21 = x2 - x1;
        double x23 = x2 - x3;
        double f1 = x21*(y2 - y3);
        double f2 = x23*(y2 - y1);
        double numerator = x21*f1 - x23*f2;
        double denominator = f1 - f2;
        return x2 - numerator/(2*denominator);
    }

    /** Find minimum along a line. */
    private ContourProperties minimise1D(Contour x0, double val, DContour dir)
    {
        final double precision1 = 1e-8;
        final double precision2 = 1e-15;

        final double goldenCut = 1.618034;
        final double initialStepSize = 0.1;
        final double maxDist = 10;

        double ya, yb, yc=0, da, db, dc;
        int step;
        ContourProperties tmpCp;

        Contour xa, xb, xc, xm;
        da = 0.0;
        xa = x0;
        ya = val;
        db = initialStepSize;
        xb = x0.plusLin(db, dir);
        tmpCp = new ContourProperties(xb);
        yb = tmpCp.getFitDistance();

        if (yb > ya) {//exchange a and b so we can assume yb is downhill
            yc = ya;
            ya = yb;
            yb = yc;
            da = initialStepSize;
            db = 0.0;
        }

        boolean calcYc = true;
        dc = db + goldenCut*(db - da);
        for (step=0; step < maxIterMin; step++) {
            if (calcYc) {
                xc = x0.plusLin(dc, dir);
                tmpCp = new ContourProperties(xc);
                yc = tmpCp.getFitDistance();
            }
            calcYc = true;
            if (yc > yb) break; // minimum between xa and xc: done

            // check for parabolic interpolation minimum
            double dba = db - da;
            double dbc = db - dc;
            double fa = dba*(yb - yc);
            double fb = dbc*(yb - ya);
            double numerator = dba*fa - dbc*fb;
            double denominator = fa - fb;
            double dm;
            if (denominator/(dbc*dba*(dba-dbc)) > 0) {// minimum ?
                dm = db - numerator/(2*denominator);
            } else {
                dm = dc + goldenCut*(dc - db);
            }

            xm = x0.plusLin(dm, dir);
            tmpCp = new ContourProperties(xm);
            double ym = tmpCp.getFitDistance();
            double dlimit = db + maxDist*(dc - db);
            if ((dm - db)*(dc - dm) > 0) {// dm is between db and dc
                if (ym < yc) {// [b,m,c]
                    da = db;
                    ya = yb;
                    db = dm;
                    yb = ym;
                    break;
                } else if (ym > yb) {// [a,b,m]
                    dc = dm;
                    yc = ym;
                    break;
                }
                // no minimum, have to extend search to greater values
                dm = dc + goldenCut*(dc - db);
            } else if ((dc - dm)*(dm - dlimit) > 0) {//minimum not too far
                if (ym > yc) {// [b,c,m]
                    da = db;
                    ya = yb;
                    db = dc;
                    yb = yc;
                    dc = dm;
                    yc = ym;
                    break;
                }
                calcYc = false;
            } else if ((dm - dlimit)*(dlimit - dc) >= 0) {// beyond limit
                dm = dlimit;
            } else {
                dm = dc + goldenCut*(dc -db);
            }
            da = db;
            ya = yb;
            db = dc;
            yb = yc;
            dc = dm;
            yc = ym;
        }
        if (step >= maxIterMin)
            log.error("minimum bracketing failed");
        //else
        //    System.err.println("minimum bracketing: ["+da+", "+dc+"]");

        // reverse intervall if da > dc
        if (da > dc) {
            double dummy = da;
            da = dc;
            dc = dummy;
        }

        // localise minimum (combined method of Brent)
        final double goldC = 0.391966;
        double x = da + goldC*(dc - da);// initialize at golden section
        double e = 0.0, d = 0.0;
        double v = x;
        double w = x;
        double u;
        xm = x0.plusLin(x, dir);
        tmpCp = new ContourProperties(xm);
        double fx = tmpCp.getFitDistance();
        double fv = fx;
        double fw = fx;
        for (step=0; step < maxIterMin; step++) {
            double dm = 0.5*(da + dc);
            double tol = precision1*Math.abs(x) + precision2;
            double twotol = 2*tol;
            if (Math.abs(x-dm) < twotol - 0.5*(dc - da)) break;
            boolean parabolaStep = false;
            if (Math.abs(e) > tol) {// fit parabola
                double r = (x - w)*(fx - fv);
                double q = (x - v)*(fx - fw);
                double p = (x - v)*q - (x - w)*r;
                q = 2*(q - r);
                if (q > 0)
                    p = -p;
                else
                    q = -q;
                r = e;
                e = d;
                //System.err.println("Parab: "+(p/q));
                if (Math.abs(p) < Math.abs(0.5*q*r) &&
                    p > q*(da - x) && p < q*(dc - x)) {
                    // use result of parabolic fit
                    d = p/q;
                    u = x + d;
                    if ((u - da) < twotol || (dc - u) < twotol) {
                        // make sure u is not too close to da or dc
                        d = -tol;
                        if (x < dm) d = tol;
                    }
                    parabolaStep = true;
                }
            }
            if (!parabolaStep) {// perform golden section step
                //System.err.println("Golden section");
                if (x < dm)
                    e = dc - x;
                else
                    e = da - x;
                d = goldC*e;
            }
            // determine point u where function is to be computed
            // making sure it is not too close to x
            if (Math.abs(d) >= tol)
                u = x + d;
            else if (d > 0)
                u = x + tol;
            else
                u = x - tol;
            xm = x0.plusLin(u, dir);
            tmpCp = new ContourProperties(xm);
            double fu = tmpCp.getFitDistance();
            // update da, dc, v, w, x
            if (fu <= fx) {
                if (u < x)
                    dc = x;
                else
                    da = x;
                v = w;
                fv = fw;
                w = x;
                fw = fx;
                x = u;
                fx = fu;
            } else {
                if (u < x)
                    da = u;
                else
                    dc = u;
                if (fu <= fw || w == x) {
                    v = w;
                    fv = fw;
                    w = u;
                    fw = fu;
                } else if (fu <= fv || v == x || v == w) {
                    v = u;
                    fv = fu;
                }
            }
        }// for step
        if (step >= maxIterMin)
            log.error("1D minimum localisation failed");
        //else
        //    log.info("minimum at " + (x0.plusLin(x,dir)));

        xm = x0.plusLin(x, dir);
        return new ContourProperties(xm);
    }


    // -- Main method --

    /** Tests our command. */
    public static void main(final String... args) throws Exception {
        final String testImagePath = "/home/adrian/Programmes/plugins_ImageJ_src/Traitement_Gouttes/src/test/resources/eauContrasteMaxStack.tif";

        // Launch ImageJ as usual.
        //final ImageJ ij = net.imagej.Main.launch(args);
        final ImageJ ij = new ImageJ();
        ij.ui().showUI();

        // Open test image.
        final ServiceHelper sh = new ServiceHelper(ij.getContext());
        final IOService io = sh.loadService(DefaultIOService.class);
        final Dataset dataset = (Dataset) io.open(testImagePath);

        // create a display for the dataset
        final ImageDisplay imageDisplay =
            (ImageDisplay) ij.display().createDisplay(dataset);

        // create a rectangle
        final RectangleOverlay rectangle = new RectangleOverlay(ij.getContext());
        rectangle.setOrigin(110, 0);
        rectangle.setOrigin(60, 1);
        rectangle.setExtent(340, 0);
        rectangle.setExtent(420, 1);
        rectangle.setLineColor(Colors.HONEYDEW);
        rectangle.setLineWidth(1);

        // add the overlays to the display
        final List<Overlay> overlays = new ArrayList<Overlay>();
        overlays.add(rectangle);
        ij.overlay().addOverlays(imageDisplay, overlays);
        for (final net.imagej.display.DataView view : imageDisplay) {
            if (view instanceof net.imagej.display.OverlayView) {
                view.setSelected(true);
            }
        }

        // display the dataset
        ij.ui().show(imageDisplay);

        // Launch the "CommandWithPreview" command.
        ij.command().run(Goutte_pendante.class, true);
    }

    // -- private classes and helper methods --

    /** An object describing a real-valued polynome in one variable. */
    private class Polynome {
        double[] coeff;

        /** Construct a polynome with given coefficients (coeff[k] is
         * the coefficient of x^k). The array must be non-null. */
        Polynome(double coeff[]) {
            if (coeff == null)
                throw new IllegalArgumentException("Cannot create a polynome with null coefficients.");
            this.coeff = coeff;
        }

        /** Evaluate polynome at x.
         *
         * @return sum_k coeff[k] * x^k
         */
        double getValueAt(double x) {
            double xp = 1, y = 0;
            for (int i=0; i<coeff.length; i++) {
                y += coeff[i]*xp;
                xp *= x;
            }
            return y;
        }

        /** Get indicated coefficient. */
        double getCoeff(int i) {
            return coeff[i];
        }

        /** Return textual description of this polynome. */
        public String toString() {
            return "polynome coeffs: "+Arrays.toString(coeff);
        }
    }

    /* Calculate the least squares linear fit to the given data
     * points. Points are ignored if either value is NaN.
     *
     * @throws IllegalArgumentException if input arrays have different lengths
     */
    private Polynome linearFit(double[] x, double[] y) {
        if (x.length != y.length)
            throw new IllegalArgumentException("linearFit: input arrays of different length (" + x.length + " != " + y.length + ")");
        double xi1 = 0, xi2 = 0, zeta0 = 0, zeta1 = 0;
        int N = 0;
        for (int i=0; i<x.length; i++) {
            if (!Double.isNaN(x[i]) && !Double.isNaN(y[i])) {
                xi1 += x[i];
                xi2 += x[i]*x[i];
                zeta0 += y[i];
                zeta1 += y[i]*x[i];
                N++;
            }
        }
        double det = xi1*xi1 - N*xi2;
        if (det == 0) {
            log.error("linear fit failed because of nil determinant for" +
                      "points:\n" +
                      Arrays.toString(x) + "\n" + Arrays.toString(y));
            return null;
        }
        double[] coeff = new double[2];
        coeff[0] = (-xi2*zeta0 + xi1*zeta1) / det;
        coeff[1] = (xi1*zeta0 - N*zeta1) / det;

        return new Polynome(coeff);
    }

    /* Calculate the least squares quadratic fit to the given data
     * points. Points are ignored if either value is NaN.
     *
     * @throws IllegalArgumentException if input arrays have different lengths
     */
    private Polynome quadraticFit(double[] x, double[] y) {
        if (x.length != y.length)
            throw new IllegalArgumentException("quadraticFit: input arrays of different length (" + x.length + " != " + y.length + ")");
        double xi1 = 0, xi2 = 0, xi3 = 0, xi4 = 0;
        double zeta0 = 0, zeta1 = 0, zeta2 = 0;
        int N = 0;
        for (int i=0; i<x.length; i++) {
            if (!Double.isNaN(x[i]) && !Double.isNaN(y[i])) {
                xi1 += x[i];
                xi2 += x[i]*x[i];
                xi3 += x[i]*x[i]*x[i];
                xi4 += x[i]*x[i]*x[i]*x[i];
                zeta0 += y[i];
                zeta1 += y[i]*x[i];
                zeta2 += y[i]*x[i]*x[i];
                N++;
            }
        }
        double m11 = N*xi4 - xi2*xi2;
        double m12 = xi1*xi2 - N*xi3;
        double m22 = N*xi2 - xi1*xi1;
        double z1 = N*zeta1 - xi1*zeta0;
        double z2 = N*zeta2 - xi2*zeta0;
        double det = m11*m22 - m12*m12;
        if (det == 0) {
            log.error("quadratic fit failed because of nil subdeterminant" +
                      "for points:\n" +
                      Arrays.toString(x) + "\n" + Arrays.toString(y));
            return null;
        }
        double[] coeff = new double[3];
        coeff[2] = (m12*z1 + m22*z2) / det;
        coeff[1] = (m11*z1 + m12*z2) / det;
        coeff[0] = (zeta0 - xi1*coeff[1] - xi2*coeff[2]) / N;

        return new Polynome(coeff);
    }

    /** A structure regrouping the geometrical parameters
     * characterizing an axisymmetric drop contour.
     */
    private class Contour {
        double tip_radius;
        double capillary_length;
        double tip_x;
        double tip_y;
        double gravity_deg;

        Contour(double tip_radius, double capillary_length,
                double tip_x, double tip_y, double gravity_deg) {
            this.tip_radius = tip_radius;
            this.capillary_length = capillary_length;
            this.tip_x = tip_x;
            this.tip_y = tip_y;
            this.gravity_deg = gravity_deg;
        }

        // -- generic accessor --

        /** get the parameter through its name */
        double get(SHAPE_PARAM param) {
            switch (param) {
            case TIP_RADIUS:
                return tip_radius;
            case CAPILLARY_LENGTH:
                return capillary_length;
            case TIP_X:
                return tip_x;
            case TIP_Y:
                return tip_y;
            case GRAVITY_ANGLE:
                return gravity_deg;
            default:
                throw new IllegalArgumentException("Don't know parameter "+param);
            }
        }

        // -- explicit accessors --

        double getTipRadius() { return tip_radius; }
        double getCapillaryLength() { return capillary_length; }
        double getTipX() { return tip_x; }
        double getTipY() { return tip_y; }
        double getGravityAngle() { return gravity_deg; }

        // -- Contour arithmetic --

        /** Build new Contour from this one translated by increment. */
        Contour plus(DContour increment) {
            return plusLin(1.0, increment);
        }

        /** Build new Contour from this one translated by a constant
         * times an increment.
         *
         * @return new(this + A * increment)
         */
        Contour plusLin(double A, DContour increment) {
            return new
                Contour(tip_radius + A*increment.get(SHAPE_PARAM.TIP_RADIUS),
                        capillary_length
                        + A*increment.get(SHAPE_PARAM.CAPILLARY_LENGTH),
                        tip_x + A*increment.get(SHAPE_PARAM.TIP_X),
                        tip_y+ A*increment.get(SHAPE_PARAM.TIP_Y),
                        gravity_deg
                        + A*increment.get(SHAPE_PARAM.GRAVITY_ANGLE));
        }

        // -- Produce readable self-documentation --

        @Override
        public String toString() {
            return "Contour[tR=" + tip_radius +
                ", cl=" + capillary_length +
                ", tX=" + tip_x +
                ", tY=" + tip_y +
                ", ang=" + gravity_deg + "]";
        }
    }

    /** A structure defining a vector in contour shape parameter
     * space. In other words, this is a difference between two
     * Contour(s)
     *
     * @see Contour
     */
    private class DContour {
        double delta_tip_radius;
        double delta_capillary_length;
        double delta_tip_x;
        double delta_tip_y;
        double delta_gravity_deg;

        DContour(double delta_tip_radius, double delta_capillary_length,
                 double delta_tip_x, double delta_tip_y,
                 double delta_gravity_deg) {
            this.delta_tip_radius = delta_tip_radius;
            this.delta_capillary_length = delta_capillary_length;
            this.delta_tip_x = delta_tip_x;
            this.delta_tip_y = delta_tip_y;
            this.delta_gravity_deg = delta_gravity_deg;
        }

        /** Construct a new DContour pointing from a to b. */
        DContour(Contour a, Contour b) {
            for (SHAPE_PARAM param: SHAPE_PARAM.values())
                put(param, b.get(param) - a.get(param));
        }

        /** Get the parameter by name */
        double get(SHAPE_PARAM param) {
            switch (param) {
            case TIP_RADIUS:
                return delta_tip_radius;
            case CAPILLARY_LENGTH:
                return delta_capillary_length;
            case TIP_X:
                return delta_tip_x;
            case TIP_Y:
                return delta_tip_y;
            case GRAVITY_ANGLE:
                return delta_gravity_deg;
            default:
                throw new IllegalArgumentException("Don't know parameter "+param);
            }
        }

        /** Modify a parameter by name. */
        void put(SHAPE_PARAM param, double newValue) {
            switch (param) {
            case TIP_RADIUS:
                delta_tip_radius = newValue;
                break;
            case CAPILLARY_LENGTH:
                delta_capillary_length = newValue;
                break;
            case TIP_X:
                delta_tip_x = newValue;
                break;
            case TIP_Y:
                delta_tip_y = newValue;
                break;
            case GRAVITY_ANGLE:
                delta_gravity_deg = newValue;
                break;
            default:
                throw new IllegalArgumentException("Don't know parameter "+param);
            }
        }

        @Override
        public String toString() {
            return "DContour[DtR=" + delta_tip_radius +
                ", Dcl=" + delta_capillary_length +
                ", DtX=" + delta_tip_x +
                ", DtY=" + delta_tip_y +
                ", Dang=" + delta_gravity_deg + "]";
        }
    }


    /** A structure regrouping various calculated properties of a
     * given contour. The parameters describing the shape of this
     * contour have to be provided to the constructor.
     */
    private class ContourProperties {
        private Contour geometry;
        private ArrayList<Point2D.Double> halfAdimProfile;
        private Path2D closedAdimProfile;
        private Shape dimShape;
        private double fitDistance;
        private double volume;
        private double surface;

        ContourProperties(Contour geometry,
                          ArrayList<Point2D.Double> halfAdimProfile,
                          Path2D closedAdimProfile,
                          Shape dimShape,
                          double fitDistance,
                          double volume,
                          double surface) {
            this.geometry = geometry;
            this.halfAdimProfile = halfAdimProfile;
            this.closedAdimProfile = closedAdimProfile;
            this.dimShape = dimShape;
            this.fitDistance = fitDistance;
            this.volume = volume;
            this.surface = surface;
        }

        ContourProperties(Contour c) {
            geometry = c;
            if (c.tip_radius > 0 && c.capillary_length > 0) {
                halfAdimProfile =
                    calculateProfile(c.tip_radius / c.capillary_length,
                                     1.5 * bounds.height * pixel_size /
                                     c.capillary_length);
                closedAdimProfile = makeClosedPath(halfAdimProfile);
                dimShape = contourToScreen(closedAdimProfile,
                                           c.capillary_length / pixel_size,
                                           c.tip_x / pixel_size,
                                           c.tip_y / pixel_size,
                                           c.gravity_deg,
                                           getBoundsArea());
                fitDistance = calcFitDistance(dimShape);
                volume = -1;
                surface = -1;
            } else {
                halfAdimProfile = null;
                closedAdimProfile = null;
                dimShape = null;
                fitDistance = Double.POSITIVE_INFINITY;
                volume = 0;
                surface = 0;
            }
        }

        public Contour getContour() {
            return geometry;
        }

        public ArrayList<Point2D.Double> getHalfAdimProfile() {
            return halfAdimProfile;
        }

        public Path2D getClosedAdimProfile() {
            return closedAdimProfile;
        }

        public Shape getDimShape() {
            return dimShape;
        }

        public double getFitDistance() {
            return fitDistance;
        }

        public double getVolume() {
            if (volume < 0) {
                final double scale = geometry.capillary_length;
                final double yMax = (geometry.tip_y - bounds.y * pixel_size)
                    / scale;
                volume = scale*scale*scale * calcVolume(halfAdimProfile, yMax);
            }
            return volume;
        }

        public double getSurface() {
            if (surface < 0) {
                final double scale = geometry.capillary_length;
                final double yMax = (geometry.tip_y - bounds.y * pixel_size)
                    / scale;
                surface = scale*scale * calcSurface(halfAdimProfile, yMax);
            }
            return surface;
        }
    }
}
