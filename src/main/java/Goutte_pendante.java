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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.HashMap;

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

import org.scijava.command.Command;
import org.scijava.command.Previewable;
import org.scijava.io.DefaultIOService;
import org.scijava.io.IOService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.service.ServiceHelper;
import org.scijava.util.Colors;
import org.scijava.util.IntRect;
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

    // The 'min' attribute requires a String, but
    // Double.toString(Double.MIN_VALUE) is not a constant to the
    // compiler, so we use an explicit value close to the real
    // constant
    @Parameter(label = "tip_radius of curvature",
               persist = false, min = "1e-300")
    private double tip_radius;

    @Parameter(label = "capillary length",
               persist = false, min = "1e-300")
    private double capillary_length;

    @Parameter(label = "tip_x coordinate",
               persist = false)
    private double tip_x;

    @Parameter(label = "tip_y coordinate",
               persist = false)
    private double tip_y;

    @Parameter(label = "gravity angle (deg)",
               persist = false)
    private double gravity_deg;

    @Parameter(label = "pixel size",
               initializer = "initPixelSize",
               persist = false,
               min = "1e-300")
    private double pixel_size;

    @Parameter(label = "density contrast times g")
    private double rho_g;

    @Parameter(label = "Surface tension",
               persist = false,
               visibility = org.scijava.ItemVisibility.MESSAGE)
    private double surface_tension = 0;

    @Parameter(label = "Adjustment penalty",
               persist = false,
               visibility = org.scijava.ItemVisibility.MESSAGE)
    private double adj_penalty = 0;

    @Parameter(label = "fit parameters checked below",
               callback = "fitButtonCB")
    private Button fitButton;

    @Parameter(label = "tip radius",
               callback = "doNothing")
    private boolean fit_include_tip_radius = true;

    @Parameter(label = "capillary length",
               callback = "doNothing")
    private boolean fit_include_capillary_length = true;

    @Parameter(label = "tip x coordinate",
               callback = "doNothing")
    private boolean fit_include_tip_x = true;

    @Parameter(label = "tip y coordinate",
               callback = "doNothing")
    private boolean fit_include_tip_y = true;

    @Parameter(label = "gravity angle",
               callback = "doNothing")
    private boolean fit_include_gravity_angle = true;

    // -- Other fields --

    /** The dimensional parameters. */
    //private HashMap paramWithDim = null;

    /** A rectangle representation of dropRegion with integer precision. */
    private IntRect bounds;

    /** Arrays containing x coordinates of the drop interface. */
    private double[] leftBorder;
    private double[] rightBorder;

    // -- Constants --

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

    // -- Command methods --

    @Override
    public void run() {
        prepareFitParams();

        ImageStack stack = imp.getStack();
        // TODO: create results table

        for (int n=0; n<stack.getSize(); n++) {
            analyseImage(stack.getProcessor(n+1));
            // TODO: copy parameters into results table

        }
    }

    // -- Previewable methods --

    @Override
    public void preview() {
        updateOverlay();
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
        bounds = new IntRect((int) Math.round(r.x),
                             (int) Math.round(r.y),
                             (int) Math.round(r.width),
                             (int) Math.round(r.height));
        //log.info("drop region: +" + bounds.x + " +" +  r.y
        //         + ", " +  r.width + " x " + r.height);

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

        findDropBorders(ip);
        int tip = leftBorder.length - 1;

        // tip curvature: linear fit to half width squared near tip
        double[] radiusSquare = new double[tipNeighbourhood];
        double[] yCoord = new double[tipNeighbourhood];
        for (int i=0; i<tipNeighbourhood; i++) {
            double halfWidth = 0.5*( rightBorder[tip-i] - leftBorder[tip-i] );
            radiusSquare[i] = halfWidth * halfWidth;
            yCoord[i] = i;
        }
        Polynome tipFit = linearFit(yCoord, radiusSquare);
        tip_radius = 0.5 * Math.abs(tipFit.getCoeff(1)) * pixel_size;

        // direction of gravity: tilt of center line
        //
        // link middle of top line to estimated center of drop
        // (one radius above tip)
        final double dX = 0.5*( rightBorder[0] + leftBorder[0]
                                - rightBorder[tip] - leftBorder[tip] );
        final double dY = tip - tip_radius / pixel_size;
        final double tiltAngleRad = Math.atan2(-dX,dY);
        gravity_deg = tiltAngleRad*180/Math.PI;

        // drop tip
        tip_x = bounds.x +
            0.5*( rightBorder[tip] + leftBorder[tip] ) * pixel_size;
        tip_y = bounds.y +
            tip + tipFit.getCoeff(0) / tipFit.getCoeff(1) * pixel_size;
        // the former expression assumes the drop is vertical, but we can
        // correct for tilt:
        tip_x -= tip_radius * Math.sin(tiltAngleRad);
        tip_y -= tip_radius * (1 - Math.cos(tiltAngleRad));

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
    }

    /** Detect drop borders and store positions in the
     * left/rightBorder arrays.
     */
    private void findDropBorders(ImageProcessor ip) {
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
        //log.info("fitButtonCB !");
        log.info("will fit:");
        if (fit_include_tip_radius) log.info("  tip radius");
        if (fit_include_capillary_length) log.info("  capillary_length");
        if (fit_include_tip_x) log.info("  tip x coordinate");
        if (fit_include_tip_y) log.info("  tip y coordinate");
        if (fit_include_gravity_angle) log.info("  gravity angle");
        prepareFitParams();
        analyseImage(imp.getProcessor());
        updateOverlay();
    }

    /** This callback is used on parameters whose change need not
     * trigger a preview() call.
     */
    private void doNothing() {
        log.info("quiet please!");
    }

    // -- Processing --

    private void prepareFitParams() {
    }

    private void updateOverlay() {
        log.info("updating overlay");
    }

    public void analyseImage(ImageProcessor ip) {
        log.info("processing ip: "+ip.toString());
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

    private class Polynome {
        double[] coeff;

        Polynome(double coeff[]) {
            this.coeff = coeff;
        }

        double getValueAt(double x) {
            double xp = 1, y = 0;
            for (int i=0; i<coeff.length; i++) {
                y += coeff[i]*xp;
                xp *= x;
            }
            return y;
        }

        double getCoeff(int i) {
            return coeff[i];
        }

        public String toString() {
            return "polynome coeffs: "+Arrays.toString(coeff);
        }
    }

    /* Calculate the least squares linear fit to the given data
     * points. Points are ignored if either value is NaN.
     *
     * @throws IllegalArgumentException if input arrays have different lengths
     */
    Polynome linearFit(double[] x, double[] y) {
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
    Polynome quadraticFit(double[] x, double[] y) {
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
}
