/**
 * Pendant drop profile integration and fitting
 *
 * Description: The pendant drop method is commonly used to measure
 * surface tensions of liquids. It consists in analysing the shape of
 * a drop hanging typically from a capillary tube and about to detach
 * (sometimes the inverse situation of a bubble forming at the bottom
 * of a liquid is preferred). The shape is very sensitive to the
 * unknown interfacial tension. The drop profile is described by only
 * one non-dimensional parameter (tip radius over capillary length),
 * although in practice five dimensional parameters can be adjusted
 * within this plug-in: tip position and curvature, tilt of symetry
 * axis and capillary length. The surface tension is calculated from
 * the latter if the density difference is given.
 *
 * This plug-in allows interactive adjustment of a profile to an image
 * of a pendant drop. An estimate of the quality of the fit is logged
 * to ImageJ's log window. The plug-in can also be asked to improve
 * the fit by varying one or several of the parameters automatically
 * (currently by successive minimisation along chosen directions
 * according to Powell).
 *
 * In principle this plug-in can also be used to estimate the surface
 * tension from the shape of a sessile drop.
 *
 * Prerequisites:
 *
 * The plug-in is run on a high contrast image of a pendant drop.
 * Pixel values are expected to be close to zero within the drop, and
 * close to saturation (i.e. 255 for 8-bit images) outside. This is
 * typically obtained by taking the picture of the drop in front of a
 * light background far away from the drop. Bright spots well within
 * the drop are not a problem, but contrast should be good in the
 * vicinity of the contour. Avoid however thresholded (binary) black
 * and white images, which generally cause the minimisation algorithm
 * to perform poorly. A little bit of blurring on such images will
 * improve the fitting performance.
 *
 * Draw a rectangular Region Of Interest (ROI) around the pendant drop
 * before calling the plug-in. The ROI should not include the inlet
 * tube, but only the free surface of the drop. A few ten or so pixels
 * margin around the drop should be sufficient. Including more of the
 * outside region will not do harm, but slow down calculations which
 * take a time proportional to the area of the ROI.
 *
 * Using the plug-in:
 *
 * Fill in any parameter values you know, and check the plausibility
 * of the initial guesses for the others. If the image was calibrated
 * for spatial scale, the corresponding scale is proposed by the
 * plug-in (but beware that ImageJ's dialogs round values to the
 * number of digits shown, so double-check if your pixel size is much
 * smaller than unity) and all length scales are expressed in
 * calibrated units. Checking the preview box will cause the profile
 * to be shown as an overlay.
 *
 * a) modify the parameters so as to improve the visual accordance of
 * theoretical profile and image. The adjustment penalty value shown
 * serves as a quantitative measure of the agreement: it should be
 * made as small as possible.
 *
 * b) check parameters to adjust and click on 'fit' for automatic
 * minimisation. The overlay is updated accordingly. Note that the
 * minimisation algorithm can get trapped in a local minimum, so try
 * starting from different initial values and check whether the
 * resulting profile looks satisfactory.
 *
 * Exit the plug-in by clicking on either OK or Cancel.
 *
 * [See more detailed description in included PDF documentation]
 */

import ij.ImagePlus;
import ij.IJ;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;
import ij.process.ByteProcessor;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.gui.GenericDialog;
import ij.gui.DialogListener;
import ij.measure.Calibration;
import java.awt.Color;
import java.awt.Panel;
import java.awt.Button;
import java.awt.Shape;
import java.awt.Rectangle;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.awt.geom.Point2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Area;
//import java.awt.geom.GeneralPath;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.TextField;
import java.awt.BorderLayout;
import java.awt.Checkbox;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;

import java.util.Vector;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.HashMap;

import java.net.URI;
import java.net.URL;
import java.net.URISyntaxException;
import java.io.InputStream;
import java.io.FileOutputStream;
import java.io.File;
import java.io.IOException;

import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.SwingConstants;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.BorderFactory;


public class Goutte_pendante implements ExtendedPlugInFilter,
                                        ActionListener, DialogListener {

    /* dialog parameters and fit parameters:
     * dialogParam[0] ... capillary length
     * dialogParam[1] ... tip x coordinate
     * dialogParam[2] ... tip y coordinate
     * dialogParam[3] ... tip radius of curvature
     * dialogParam[4] ... deviation of gravity from vertical in degrees
     * dialogParam[5] ... density contrast times gravitational acceleration (rho*g)
     * dialogParam[6] ... length scale (pixel size)
     *
     * all lengths in dialogParam[] are in calibrated units, if any,
     * otherwise in pixels (in which case dialogParam[6] == 1)
     *
     * fitParam[0] ... capillary length
     * fitParam[1] ... tip x coordinate
     * fitParam[2] ... tip y coordinate
     * fitParam[3] ... tip radius curvature
     * fitParam[4] ... gravity angle plus PI (radians)
     *
     * all lengths in fitParam[] have the same units as in dialogParam[];
     * the angle is in radians and has 2*PI added to be O(1).
     */
    private final static String[] dialogNames = { "capillary length",
                                    "tip_x coordinate",
                                    "tip_y coordinate",
                                    "tip_radius of curvature",
                                    "gravity angle (deg)",
                                    "density contrast times g",
                                    "scale (pixel size)" };
    private final static int NdialogParam = dialogNames.length;
    //@GuardedBy dialog
    private double[] dialogParam = new double[NdialogParam];
    private final static String[] fitNames = { "capillary length",
                                       "tip_x coordinate",
                                       "tip_y coordinate",
                                       "tip_radius of curvature",
                                       "gravity angle" };
    private final static int NfitParam = fitNames.length;
    private boolean[] fitMe = new boolean[NfitParam];
    
    private final static String pluginMenuName = "Pendant drop";

    private ImagePlus imp;

    private GenericDialog dialog;
    private Button fitButton;
    private final static String surfaceTensionString = "Surface tension: ";
    private java.awt.Label surfaceTensionLabel;
    private final static String adjustmentPenaltyString = "Fit penalty: ";
    private java.awt.Label adjustmentPenaltyLabel;

    // for the 'about' box
    private final java.util.Map<JButton,URI> uris = 
        new java.util.LinkedHashMap<JButton,URI>();

    
    /** Plug-in starts here, and ends here. At the beginning, just
     * does basic tests. At the end, prints results into a table.
     *
     * @see ij.plugin.ExtendedPlugInFilter interface description.
     */
    public int setup(String arg,
                     ImagePlus imp) {

        if (arg.equals("about")) {
            showAbout();
            return DONE;
        } else if (arg.equals("final")) {
            IJ.log("__Drop shape summary__");
            IJ.log("capillary length: "+dialogParam[0]);
            IJ.log("tip x coordinate: "+dialogParam[1]);
            IJ.log("tip y coordinate: "+dialogParam[2]);
            IJ.log("tip radius of curvature: "+dialogParam[3]);
            IJ.log("deviation gravity from vertical (deg): "+dialogParam[4]);
            IJ.log("density contrast * grav. accel. (rho*g): "+dialogParam[5]);
            IJ.log("pixel size: "+dialogParam[6]+"\n");
            IJ.log("surface tension: "+dialogParam[0]*dialogParam[0]*dialogParam[5]);
        }

        if (IJ.versionLessThan("1.43u")) {
            IJ.error(pluginMenuName,
                     "This plug-in requires ImageJ version >= 1.43u.");
            return DONE; 
        }
        if (imp == null) {
            IJ.noImage();
            return DONE;
        }
        this.imp = imp;

        return DOES_8G | DOES_16 | DOES_32 | ROI_REQUIRED |
            PARALLELIZE_STACKS | FINAL_PROCESSING;
    }

    /** Called after setup().
     *
     * @see ij.plugin.ExtendedPlugInFilter interface description.
     */
    public int showDialog(ImagePlus imp,
                          String command,
                          PlugInFilterRunner pfr) {

        // Estimate reasonable parameters from the current image
        estimateParameters(imp.getProcessor());

        dialog = new GenericDialog(command);
	setDialogActive(true);
        dialog.addDialogListener(this);
        for (int n=0; n < NdialogParam; n++) {
            dialog.addNumericField(dialogNames[n], dialogParam[n], 8);
        }
        dialog.addMessage(surfaceTensionString + "<no value>");
        surfaceTensionLabel = (java.awt.Label)dialog.getMessage();
        dialog.addMessage(adjustmentPenaltyString + "<no value>");
        adjustmentPenaltyLabel = (java.awt.Label)dialog.getMessage();
        dialog.addMessage("Fit the following parameters:");
        dialog.addCheckboxGroup(NfitParam, 1, fitNames, fitMe);
        Panel panel = new Panel();
        fitButton = new Button("Fit now");
        fitButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    fitButtonPressed();
                }
            });
        panel.add(fitButton);
        dialog.addPanel(panel);
        dialog.addPreviewCheckbox(pfr);
        dialog.getPreviewCheckbox().setState(true);

        dialog.showDialog();
	setDialogActive(false);

        if (dialog.wasCanceled()) return DONE;

        return DOES_8G | DOES_16 | DOES_32 | ROI_REQUIRED |
            PARALLELIZE_STACKS | FINAL_PROCESSING;
    }

    /** run() should behave differently according to whether it is run
     * to obtain a preview, or for the final image, so we provide a
     * flag that tells the context of the call
     */
    private boolean dialogActive;
    synchronized void setDialogActive(boolean status) {
        dialogActive = status;
    }
    synchronized boolean isDialogActive() {
        return dialogActive;
    }

    /** Checks the validity of the parameters every time the dialog is
     * modified, and copies them to the dialogParam array. The return
     * value determines whether a potential preview calculation may be
     * triggered or not.
     *
     * @param e The event that has been generated by the user action
     * in the dialog. Note that e is null if the dialogItemChanged
     * method is called after the user has pressed the OK button or if
     * the GenericDialog has read its parameters from a macro.
     *
     * @return True if the dialog input is valid. False
     * disables the OK button and preview (if any).
     *
     * @see ij.gui.DialogListener interface
     */
    public boolean dialogItemChanged(GenericDialog gd, java.awt.AWTEvent e) {
        boolean change = false;
        boolean haveFreeParam = false;

        synchronized (dialog) {
            for (int n=0; n < NdialogParam; n++) {
                double val = gd.getNextNumber();
                if (Double.isNaN(val)) return false;
                // lengths and scale must be positive
                if ((n < 4 || n == 6) && val <= 0) return false;
                if (val != dialogParam[n]) change |= true;
                dialogParam[n] = val;
            }

            for (int n=0; n < NfitParam; n++) {
                fitMe[n] = gd.getNextBoolean();
                haveFreeParam |= fitMe[n];
            }
   
        }
        if (!haveFreeParam) {
            fitButton.setEnabled(false);
            return false;
        } else {
            fitButton.setEnabled(true);
        }
        
        return change;
    }


    /** Set fit parameters from dialog values. */
    void setDropFromDialog(DropInfo di) {
        double[] p = di.getFitParam();
        synchronized (dialog) {
            // set fit parameters from dialog params
            p[0] = dialogParam[0];
            p[1] = dialogParam[1];
            p[2] = dialogParam[2];
            p[3] = dialogParam[3];
            p[4] = (dialogParam[4]/180+2)*Math.PI;
            di.setDensityGrav(dialogParam[5]);
            di.setPixelScale(dialogParam[6]);
        }
    }
    
    
    /** Set dialog parameters from fitted values. */
    @SuppressWarnings("unchecked")
    void setDialogFromDrop(DropInfo di) {
        double[] p = di.getFitParam();
        synchronized (dialog) {
            // set dialog parameters from fit parameters
            dialogParam[0] = p[0];
            dialogParam[1] = p[1];
            dialogParam[2] = p[2];
            dialogParam[3] = p[3];
            dialogParam[4] = ((p[4]/Math.PI-2)*180);
            //dialogParam[5] = di.getDensityGrav();
            //dialogParam[6] = di.getPixelScale();

            Vector<TextField> nf = (Vector<TextField>)dialog.getNumericFields();
            for (int n=0; n < NdialogParam; n++) {
                nf.elementAt(n).setText(Double.toString(dialogParam[n]));
            }

            adjustmentPenaltyLabel.setText(adjustmentPenaltyString + di.getQuality());
            surfaceTensionLabel.setText(surfaceTensionString +
                                        p[0]*p[0]*dialogParam[5]);

        }
    }

    
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
            throw new IllegalArgumentException("linearFit: input arrays of different length ("+x.length+" != "+y.length+")");
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
            IJ.log("linear fit failed because of nil determinant for points:\n"
                   +Arrays.toString(x)+"\n"+Arrays.toString(y));
            return null;
        }
        double[] coeff = new double[2];
        coeff[0] = (-xi2*zeta0 + xi1*zeta1)/det;
        coeff[1] = (xi1*zeta0 - N*zeta1)/det;

        return new Polynome(coeff);
    }

    /* Calculate the least squares quadratic fit to the given data
     * points. Points are ignored if either value is NaN.
     *
     * @throws IllegalArgumentException if input arrays have different lengths
     */
    Polynome quadraticFit(double[] x, double[] y) {
        if (x.length != y.length)
            throw new IllegalArgumentException("quadraticFit: input arrays of different length ("+x.length+" != "+y.length+")");
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
            IJ.log("linear fit failed because of nil subdeterminant for points:\n"
                   +Arrays.toString(x)+"\n"+Arrays.toString(y));
            return null;
        }
        double[] coeff = new double[3];
        coeff[2] = (double)((m12*z1 + m22*z2)/det);
        coeff[1] = (double)((m11*z1 + m12*z2)/det);
        coeff[0] = (double)((zeta0 - xi1*coeff[1] - xi2*coeff[2])/N);

        return new Polynome(coeff);
    }

    /** Try to guess sensible initial values for the parameters. */
    void estimateParameters(ImageProcessor ip) {
        // pixel size: get from image if set
        Calibration cal = imp.getCalibration();
        if (cal.scaled())
            dialogParam[6] = (double)cal.pixelWidth;
        else
            dialogParam[6] = 1;

        // gather rough info on the drop
        Rectangle roi = ip.getRoi().getBounds();
        final double threshold = 128.f;
        int tip=0;
        double[] centerX = new double[roi.height];
        double[] centerY = new double[roi.height];
        double[] halfWidth = new double[roi.height];
        for (int y=0; y < roi.height; y++) {// y increases upwards in ROI
            centerY[y] = roi.y + roi.height - 1 - y;
            double centerXacc = 0;// accumulator for center of gravity along line
            int count = 0;
            for (int x=0; x<roi.width; x++) {
                if (ip.getPixelValue(roi.x + x,
                                     roi.y + roi.height - 1 - y) < threshold) {
                    count++;
                    centerXacc += x;
                }
            }
            if (count > 0) {
                centerX[y] = roi.x + centerXacc/count;
                if (tip == 0)
                    tip = y;
                halfWidth[y] = 0.5f*count;
            } else {
                centerX[y] = Double.NaN;
                halfWidth[y] = Double.NaN;
            }
        }

        // tip curvature: linear fit to halfWidth near tip
        final int tipNeighbourhood = 5;// fit on these many points
        double[] radiusSquare = new double[tipNeighbourhood];
        for (int i=0; i<tipNeighbourhood; i++) {
            radiusSquare[i] = halfWidth[tip+i]*halfWidth[tip+i];
        }
        Polynome tipFit = linearFit(Arrays.copyOfRange(centerY, tip,
                                                       tip+tipNeighbourhood),
                                    radiusSquare);
        dialogParam[3] = Math.abs(tipFit.getCoeff(1))/2*dialogParam[6];

        // direction of gravity: tilt of center line
        //
        // link middle of top line to estimated center of drop
        // (one radius above tip)
        final double axisX = centerX[roi.height-1] - centerX[tip];
        final double axisY = centerY[roi.height-1] - centerY[tip]
            + dialogParam[3]/dialogParam[6];
        final double tiltAngleRad = Math.atan2(axisX,-axisY);
        dialogParam[4] = (double)(tiltAngleRad*180/Math.PI);

        // detected drop tip (more precise than centerY[tip]*dialogParam[6])
        dialogParam[1] = centerX[tip]*dialogParam[6];
        dialogParam[2] = -tipFit.getCoeff(0)/tipFit.getCoeff(1)*dialogParam[6];
        // the former expression assumes the drop is vertical, but we can
        // correct for tilted image:
        dialogParam[1] -= dialogParam[3] * (double)Math.sin(tiltAngleRad);
        dialogParam[2] -= dialogParam[3] * (1 - (double)Math.cos(tiltAngleRad));

        // capillary length: from curvature difference

        final int bellyNeighbourhood = (roi.height-tip)/8;// fit on up to twice
        // these many points

        // find maximum thickness
        int yBelly = tip;
        for (int y=tip+1; y < roi.height; y++) {
            if (!Double.isNaN(halfWidth[y]))
                if (halfWidth[y] > halfWidth[yBelly])
                    yBelly = y;
        }
        if (roi.height - yBelly < bellyNeighbourhood)// too close to top;
            yBelly = roi.height - 1 - bellyNeighbourhood;

        //IJ.log("tip = "+tip+", yBelly = "+yBelly+
        //       ", length = "+centerY.length);

        // fit parabola to neighbourhood of yBelly
        int yBelow = Math.max(tip, yBelly - bellyNeighbourhood);
        int yAbove = Math.min(centerY.length - 1, yBelly + bellyNeighbourhood);
        double[] xB = Arrays.copyOfRange(centerY, yBelow, yAbove);
        double[] yB = Arrays.copyOfRange(halfWidth, yBelow, yAbove);
        double mxB = 0, myB = 0;
        int cx = 0, cy = 0;
        for (int y=0; y < xB.length; y++) {
            if (!Double.isNaN(xB[y])) {
                mxB += xB[y];
                cx ++;
            }
            if (!Double.isNaN(yB[y])) {
                myB += yB[y];
                cy ++;
            }
        }
        mxB /= cx;
        myB /= cy;
        for (int y=0; y < xB.length; y++) {
            xB[y] -= mxB;
            yB[y] -= myB;
        }

        Polynome bellyFit = quadraticFit(xB, yB);

        // and calculate curvatures at yBelly
        double c1 = 1/halfWidth[yBelly];
        double c2 = - 2*bellyFit.getCoeff(2) /
            (double)(Math.cos(Math.atan(bellyFit.getCoeff(1))));
        double c0 = dialogParam[6]/dialogParam[3];// tip curvature

        //IJ.log("c0 = "+c0+", c1 = "+c1+", c2 = "+c2);

        // pressure difference between yBelly and tip gives capillary
        // length estimate
        dialogParam[0] = (double)Math.sqrt((tip-yBelly)/(c1+c2-2*c0))*dialogParam[6];

        // initialise density contrast if image calibrated
        // (not very important...)
        // for water at 4 deg C on earth etc, in grams per mm^2 and seconds^2:
        if (dialogParam[6] == 1)// length scale uncalibrated?
            dialogParam[5] = 1.0f;// don't bother
        else
            dialogParam[5] = 9.81f;
    }

    /** Registered as an ActionListener on the 'fit' button,
     * is therefore called when the latter is pressed.
     */
    void fitButtonPressed() {
        final DropInfo di = DropInfo.getDropInfo(imp.getProcessor());
        setDropFromDialog(di);
        
        runFit(di);

        showCurve(di.getShape());
        setDialogFromDrop(di);
    }

    /** Calculate new contour. Called when preview is requested and
     * any parameter is modified, or for the final processing.
     */
    public void run(ImageProcessor ip) {
        final DropInfo di = DropInfo.getDropInfo(ip);

        // decide whether to do full fit or merely update contour
        if (isDialogActive()) {
            // prepare with latest parameters from dialog
            setDropFromDialog(di);

            runContour(di);
            
            showCurve(di.getShape());
            setDialogFromDrop(di);
        } else {
            runFit(di);
        }
    }

    /** Updates the contour. */
    public void runContour(DropInfo di) {
        // calculate the drop profile
        double[] param = di.getFitParam();
        double pixelScale = di.getPixelScale();
        
        Path2D p = calculateProfile(param[3]/param[0],
                                    di.bounds.height*pixelScale/param[0]);
        Shape s = dropToScreen(p, param[0]/pixelScale, param[1]/pixelScale,
                               param[2]/pixelScale, param[4], di);

        di.setProfile(p);
        di.setShape(s);
    }

    /** Find a drop shape that minimises the adjustment penalty
     * function, by varying a subset of the shape parameters. */
    public boolean runFit(DropInfo di) {
            // check we have at least one parameter to fit
            boolean haveFreeParam = false;
            for (int i=0; i < fitMe.length; i++)
                haveFreeParam |= fitMe[i];
            if (!haveFreeParam) {
                IJ.error(pluginMenuName,
                         "Please check at least one parameter to fit");
                return false;
            }
            
            // and start the minimisation
            double[] fitParam = di.getFitParam();
            minPowell(fitParam, fitMe, calcMatchVal, 100);

            runContour(di);

            return true;
    }

    double calcDropSurface(Path2D drop) {
        final double flatness = 1e-5f; // should not have any effect
                                        // as drop is a polygon
        double oz, or, z, r; // initial and final points of segments
        final double[] coords = new double[6]; // segment coordinates
        double surface = 0;

        PathIterator iter = drop.getPathIterator(null, flatness);
        if (iter.isDone()) return 0; // should not happen: no drop points
        iter.currentSegment(coords);
        or = coords[0];
        oz = coords[1];
        iter.next();
        while (!iter.isDone()) {
            int type = iter.currentSegment(coords);
            r = coords[0];
            z = coords[1];
            if (r<0 || type == PathIterator.SEG_CLOSE) {
                break;
            }
            surface += Math.PI*(r+or)*Math.sqrt((r-or)*(r-or)+(z-oz)*(z-oz));
            or = r;
            oz = z;
            iter.next();
        }
        return surface;
    }

    /** Set given curve as overlay on image. */
    void showCurve(Shape c) {
        if (c == null) return;
        imp.setOverlay(c, Color.red, null);
    }

    /** Method of the ExtendedPlugInFilter which we do not use here.
     *
     * @see ij.plugin.ExtendedPlugInFilter interface description.
     */
    public void setNPasses(int nPasses) {}


    /** Hydrostatic pressure equilibrium equations for an axisymetric drop.
     * variables (r, psi, kappa, z) represent
     *
     * r ... drop radius in cylindrical coordinates
     * psi ... interface tilt to horizontal
     * kappa ... radial curvature
     * z ... height
     */
    private void deriv(double etat[], double deriv[]) {
        final double c = Math.cos(etat[1]);
        final double s = Math.sin(etat[1]);
        deriv[0] = c;
        deriv[1] = etat[2];
        deriv[2] = - s + c*(s/etat[0] - etat[2])/etat[0];
        deriv[3] = s;
    }

    /** Calculate non-dimensional drop profile through integration of
     * the hydrostatic pressure equilibrium. Assumes axisymetry.
     * Expects tip radius and maximum z value to be given normalised
     * by capillary length. */
    Path2D calculateProfile(double tipRadius, double zMax) {
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
        double[] etat, tmp, k1, k2, k3, k4;
        int nVar = 4;// r, psi, kappa, z
        etat = new double[nVar];
        tmp = new double[nVar];
        k1 = new double[nVar];
        k2 = new double[nVar];
        k3 = new double[nVar];
        k4 = new double[nVar];
        // intial conditions using approximate solution
        etat[0] = tipRadius*Math.sin(s/tipRadius) +
            Math.pow(s,5)/(40*tipRadius*tipRadius);// r
        etat[1] = s*(1 - 0.125f*s*s)/tipRadius;// psi
        etat[2] = (1 - 0.375f*s*s)/tipRadius;// kappa
        etat[3] = 0.5f*s*s*(1 - s*s*(0.75f + 1/(tipRadius*tipRadius))/12) /
            tipRadius;// z
        while (etat[3] < zMax) {// integrate until tall enough
            deriv(etat,k1);
            for (int i=0; i<nVar; i++) tmp[i] = etat[i] + 0.5f*ds*k1[i];
            deriv(tmp,k2);
            for (int i=0; i<nVar; i++) tmp[i] = etat[i] + 0.5f*ds*k2[i];
            deriv(tmp,k3);
            for (int i=0; i<nVar; i++) tmp[i] = etat[i] + ds*k3[i];
            deriv(tmp,k4);
            for (int i=0; i<nVar; i++)
                etat[i] = etat[i] + ds*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
            if (k4[3] < 0) break; // abort if profile bends back downwards
            profile.add(new Point2D.Double((double)etat[0],(double)etat[3]));
        }
        //System.err.println("r="+etat[0]);
        //System.err.println("z="+etat[3]);
        //System.err.println("psi="+etat[1]);
        //System.err.println("kappa="+etat[2]);

        // finally construct closed drop contour from profile by
        // concatenating the profile with its reflection about the
        // y-axis
        Path2D.Double drop = new Path2D.Double();
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
     *  @param angle deviation of gravity from vertical in degrees
     */
    Shape dropToScreen(Path2D drop, double scale, double tipX, double tipY,
                       double angle, DropInfo di) {
        //System.err.println("scale = "+scale+", tipX = "+tipX+", tipY = "+tipY+", angle = "+angle);

        // read this bottom-up to get the order of transformations right...
        AffineTransform toScreenT;
        toScreenT = AffineTransform.getTranslateInstance(tipX, tipY);
        toScreenT.scale(scale, scale);
        toScreenT.rotate(angle);
        toScreenT.quadrantRotate(2);
        //return drop.createTransformedShape(toScreenT);

        // Generate a low resolution version of the drop path. If the
        // path is polygonal with roughly pixel precision we will not
        // see the difference in the fit (except on the computation
        // time...)
        final double flatness = 0.5; // probably no effect unless the
        // flattening path iterator also
        // flattens polygons
        final double minDist2 = 1.0f; // square of minimal distance
        // between points in the low res
        // version (in pixels)
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
        visibleDrop.intersect(di.getBoundsArea());

        //System.err.println("original: "+getShapeLength(drop));
        //System.err.println("low res : "+getShapeLength(lowresDrop));
        //System.err.println("visible: "+getShapeLength(visibleDrop));

        return visibleDrop;
    }

    // used only for debugging purposes
    int getShapeLength(Shape s) {
        PathIterator iter = s.getPathIterator(null);
        int size = 0;
        while (!iter.isDone()) {
            size++;
            iter.next();
        }
        return size;
    }

    double fitQuality(Shape drop) {
        double[] leftIntersect = new double[roi.height];
        double[] rightIntersect = new double[roi.height];
        for (int y = 0; y < roi.height; y++) {
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
        //IJ.log("first segment type " + segType);
        //IJ.log("CLOSE == " + PathIterator.SEG_CLOSE);
        //IJ.log("MOVETO == " + PathIterator.SEG_MOVETO);
        //IJ.log("LINETO == " + PathIterator.SEG_LINETO);
        //IJ.log("CUBICTO == " + PathIterator.SEG_CUBICTO);
        //IJ.log("QUADTO == " + PathIterator.SEG_QUADTO);
        xo = coord[0] - roi.x;
        yo = coord[1] - roi.y;
        yof = (int)Math.floor(yo);
        //IJ.log("first point (" + xo + ", " + yo + "), roi.y = " + roi.y);
        leftIntersect[(int)yo] = xo;

        pi.next();

        while (!pi.isDone()) {
            segType = pi.currentSegment(coord);
            if (segType != PathIterator.SEG_LINETO) {
                pi.next();
                continue;
            }
            xn = coord[0] - roi.x;
            yn = coord[1] - roi.y;
            // find out if this segment crossed a horizontal
            //line of integer coordinate
            ynf = (int)Math.floor(yn);
            if (ynf != yof) {// integer bound crossed !
                int ymin = Math.max(0, Math.min(ynf, yof) + 1);
                int ymax = Math.min(roi.height - 1, Math.max(ynf, yof));
                for (int yi = ymin; yi <= ymax; yi++) {
                double frac = (yi - yo) / (yn - yo);
                double xi = xo + frac * (xn - xo);
                //IJ.log("on coupe " + yi + " pour frac = "+frac+", en x = " + xi);
                if (ynf > yof) // left side
                  leftIntersect[yi] = xi;
                else
                  rightIntersect[yi] = xi;
                }
            }
            //IJ.log("lineto (" + xn + ", " + yn + ", " + ynf + ")");
            xo = xn;
            yo = yn;
            yof = ynf;
            pi.next();
        }
        rightIntersect[(int)yo] = xo;

        // add values where there is no contour
        double dummyIntersect = roi.width/2;
        for (int y=0; y < roi.height; y++) {
            if (Double.isNaN(leftIntersect[y]) || Double.isNaN(rightIntersect[y])) {
                if (Double.isNaN(leftIntersect[y]))
                    leftIntersect[y] = dummyIntersect;
                if (Double.isNaN(rightIntersect[y]))
                    rightIntersect[y] = dummyIntersect;
            } else {
                dummyIntersect = 0.5*(leftIntersect[y] + rightIntersect[y]);
            }
            //IJ.log("y = " + y + ": left = " + leftIntersect[y] + ", right = " + rightIntersect[y]);
        }

        // calculate penalty by summing over the distances between the measured
        // and the theoretical border position.
        double Q = 0;
        double xml = 0;
        double xtl = 0;
        double xmr = 0;
        double xtr = 0;
        for (int y = 0; y < yDropTip; y++){
          if (!Double.isNaN (leftBorder[y])){
              xml = leftBorder[y];
            if (!Double.isNaN (leftIntersect[y])){
                xtl = leftIntersect[y];
            } else {
                xtl = roi.width/2;
            }
            Q += (xml - xtl)*(xml - xtl); //alternative: Math.abs(xml - xtl);
           }

        if (!Double.isNaN(rightBorder[y])){
              xmr = rightBorder[y];
              if (!Double.isNaN(rightIntersect[y])){
                  xtr = rightIntersect[y];
              } else {
                  xtr = roi.width/2;
              }
            Q += (xmr - xtr)*(xmr - xtr); //alternative: Math.abs(xmr - xtr);
            }
        }
        return Q;
    }


    private class ParameterSpacePoint {
        double[] x;
        double val;
    }

    private interface Function {
        double function(double[] x);
    }

    /* This function calculates the quality of the fit. */
    final Function calcMatchVal = new Function() {
            public double function(double[] x) {
                Path2D p=null;
                // never try to calculate profile for negative lengths
                if (x[0] < 0 || x[3] <= 0) return Double.POSITIVE_INFINITY;
                try {
                    p = calculateProfile(x[3]/x[0],
                                         roi.height*dialogParam[6]/x[0]);
                } catch (Exception e) {
                    e.printStackTrace();
                    System.err.println("roi.height="+roi.height+
                                       ", dialogParam[6]="+dialogParam[6]+", x[0]="+x[0]);
                }
                Shape s = dropToScreen(p, x[0]/dialogParam[6], x[1]/dialogParam[6],
                                       x[2]/dialogParam[6], x[4], di);
                double Q = fitQuality(s);
                //System.err.println("Q = "+Q+" at x[0] = "+x[0]+", x[1] = "+x[1]+", x[2] = "+x[2]+", x[3] = "+x[3]+", x[4] = "+x[4]);
                return Q;
            }
        };


    double minParabola(double x1, double x2, double x3,
                       double y1, double y2, double y3) {
        double x21 = x2 - x1;
        double x23 = x2 - x3;
        double f1 = x21*(y2 - y3);
        double f2 = x23*(y2 - y1);
        double numerator = x21*f1 - x23*f2;
        double denominator = f1 - f2;
        //if (denominator/(x23*x21*(x21-x23)) > 0)// minimum
        return x2 - numerator/(2*denominator);
        //else
        //    return Math.min(y1,Math.min(y2,y3));
    }

    ParameterSpacePoint minimise1D(double[] x0, double val, double[] xdir,
                                   Function minMe, int maxNumStep){
        final double precision1 = 1e-8;
        final double precision2 = 1e-15;

        final double goldenCut = 1.618034;
        final double initialStepSize = 0.1;
        final double maxDist = 10;
        double ya, yb, yc=0, da, db, dc;
        int step;

        final int N = x0.length;
        double[] xa = new double[N];
        double[] xb = new double[N];
        double[] xc = new double[N];
        double[] xm = new double[N];
        da = 0.0;
        for (int i=0; i < N; i++) xa[i] = x0[i];
        ya = val;
        db = initialStepSize;
        for (int i=0; i < N; i++) xb[i] = x0[i] + db*xdir[i];
        yb = minMe.function(xb);

        if (yb > ya) {//exchange a and b so we can assume yb is downhill
            yc = ya;
            ya = yb;
            yb = yc;
            da = initialStepSize;
            db = 0.0;
        }

        boolean calcYc = true;
        dc = db + goldenCut*(db - da);
        for (step=0; step < maxNumStep; step++) {
            if (quitLengthyCalculations()) break;// user requested interrupt
            if (calcYc) {
                for (int i=0; i < N; i++) xc[i] = x0[i] + dc*xdir[i];
                yc = minMe.function(xc);
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
                //System.err.println("minimum at dm = "+dm);
            } else {
                dm = dc + goldenCut*(dc - db);
                //System.err.println("maximum at dm = "+
                //                   (db - numerator/(2*denominator)));
            }
            //double dm = minParabola(da, db, dc, ya, yb, yc);

            for (int i=0; i < N; i++) xm[i] = x0[i] + dm*xdir[i];
            double ym = minMe.function(xm);
            //System.err.println("a ="+da+"\tb ="+db+"\tc ="+dc+"\tm ="+dm);
            //System.err.println("ya="+ya+"\tyb="+yb+"\tyc="+yc+"\tym="+ym);
            double dlimit = db + maxDist*(dc - db);
            if ((dm - db)*(dc - dm) > 0) {// dm is between db and dc
                if (ym < yc) {// [b,m,c]
                    da = db;
                    ya = yb;
                    db = dm;
                    yb = ym;
                    //System.err.println("B1");
                    break;
                } else if (ym > yb) {// [a,b,m]
                    dc = dm;
                    yc = ym;
                    //System.err.println("B2");
                    break;
                }
                // no minimum, have to extend search to greater values
                dm = dc + goldenCut*(dc - db);
                //System.err.println("B3 trying [b,c,gold]");
            } else if ((dc - dm)*(dm - dlimit) > 0) {//minimum not too far
                if (ym > yc) {// [b,c,m]
                    da = db;
                    ya = yb;
                    db = dc;
                    yb = yc;
                    dc = dm;
                    yc = ym;
                    //System.err.println("B4");
                    break;
                    //} else {
                    //    db = dc;
                    //    yb = yc;
                    //    dc = dm;
                    //    yc = ym;
                    //    dm = dc + goldenCut*(dc -db);
                }
                calcYc = false;
                //System.err.println("B5 trying [b,c,m]");
            } else if ((dm - dlimit)*(dlimit - dc) >= 0) {// beyond limit
                dm = dlimit;
                //System.err.println("B6 trying [b,c,limit]");
            } else {
                dm = dc + goldenCut*(dc -db);
                //System.err.println("B7 trying [b,c,gold]");
            }
            da = db;
            ya = yb;
            db = dc;
            yb = yc;
            dc = dm;
            yc = ym;
        }
        if (step >= maxNumStep)
            IJ.log("minimum bracketing failed");
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
        for (int i=0; i < N; i++) xm[i] = x0[i] + x*xdir[i];
        double fx = minMe.function(xm);
        double fv = fx;
        double fw = fx;
        for (step=0; step < maxNumStep; step++) {
            if (quitLengthyCalculations()) break;// user requested interrupt
            //System.err.println("left = "+da+", right = "+dc);
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
            for (int i=0; i < N; i++) xm[i] = x0[i] + u*xdir[i];
            double fu = minMe.function(xm);
            // update da, dc, v, w, x
            //System.err.println("x = "+x);
            if (fu <= fx) {
                //System.err.println("u = "+u+", x = "+x);
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
                //System.err.println("x = "+x+", u = "+u);
                if (u < x)
                    da = u;
                else
                    dc = u;
                if (fu <= fw || w == x) {
                    //System.err.println("u = "+u+", x = "+x);
                    v = w;
                    fv = fw;
                    w = u;
                    fw = fu;
                } else if (fu <= fv || v == x || v == w) {
                    //System.err.println("v = "+v+", fv = "+fv);
                    v = u;
                    fv = fu;
                }
            }
        }// for step
        if (step >= maxNumStep)
            IJ.log("minimum localisation failed");
        else {
            IJ.log("minimum for delta = ["+x*xdir[0]+", "+x*xdir[1]+", "+x*xdir[2]+", "+x*xdir[3]+", "+x*xdir[4]+"], step =  "+x);
        }

        ParameterSpacePoint p = new ParameterSpacePoint();
        for (int i=0; i < N; i++) xm[i] = x0[i] + x*xdir[i];
        p.x = xm;
        p.val = fx;

        return p;
    }

    double minPowell(double[] fitParam, boolean[] fitMe, Function minMe, int maxNumStep) {
        final double precision1 = 1e-8;
        final double precision2 = 1e-15;
        final int N = fitParam.length;
        double[] xCur = Arrays.copyOf(fitParam,N);
        double[] xOld = new double[N];
        double[] xE = new double[N];
        double fCur = minMe.function(xCur);
        double fOld, fE;
        ParameterSpacePoint pMin;
        int Nreduced = 0;
        for (int i=0; i < N; i++)
            if (fitMe[i]) Nreduced++;
        // create list of directions
        ArrayList<double[]> directions = new ArrayList<double[]>(Nreduced);
        for (int i=0; i < N; i++)
            if (fitMe[i]) {
                double[] direction = new double[N];
                Arrays.fill(direction, 0.0);
                direction[i] = 1.0;
                directions.add(direction);
            }
        int step;
        for (step=0; step < maxNumStep; step++) {
            if (quitLengthyCalculations()) break;// user requested interrupt
            fOld = fCur;
            for (int i=0; i < N; i++) xOld[i] = xCur[i];
            double deltaFmax = 0;
            int maxDir = 0;// direction of steepest descent
            int numDownSteps = 0;// number of descending steps this round
            // loop over directions and find direction which yields
            // the largest function decrease
            for (int i=0; i < Nreduced; i++) {
                pMin = minimise1D(xCur, fCur, directions.get(i),
                                  minMe, maxNumStep);
                if (quitLengthyCalculations()) break;// requested interrupt
                double deltaF = fCur - pMin.val;
                IJ.log("dir "+i+": fCur = "+fCur+": deltaF = "+deltaF+", deltaFmax = "+deltaFmax);
                if (deltaF > deltaFmax) {
                    maxDir = i;
                    deltaFmax = deltaF;
                }
                if (deltaF > 0) {
                    xCur = pMin.x;
                    fCur = pMin.val;
                    numDownSteps ++;
                }
                Path2D p = calculateProfile(xCur[3]/xCur[0],
                                            roi.height*dialogParam[6]/xCur[0]);
                Shape s = dropToScreen(p, xCur[0]/dialogParam[6], xCur[1]/dialogParam[6],
                                       xCur[2]/dialogParam[6], xCur[4]);
                showCurve(s);
            }
            // test for convergence
            IJ.log("change = "+Math.abs(fCur - fOld)+", threshold = "+(precision1*Math.abs(fCur) + precision2));
            if (Math.abs(fCur - fOld) < precision1*Math.abs(fCur) + precision2)
                break;
            // construct extrapolated point and direction from xCur and xOld
            double[] newDirection = new double[N];
            for (int i=0; i < N; i++) {
                newDirection[i] = xCur[i] - xOld[i];
            }
            if (numDownSteps > 1) {
                IJ.log("testing direction ["+newDirection[0]+", "+
                       newDirection[1]+", "+newDirection[2]+", "+
                       newDirection[3]+", "+newDirection[4]+"]");
                for (int i=0; i < N; i++) xE[i] = xCur[i] + newDirection[i];
                fE = minMe.function(xE);
                IJ.log("fOld = "+fOld+", fCur = "+fCur+", fE = "+fE);
                // now there are 3 points: (xOld, fOld), (xCur, fCur)
                // and (xE, fE); check where minimum is on this line
                if (fE < fCur) {
                    double test = 2*(fOld - 2*fCur + fE) *
                        (fOld - fCur - deltaFmax) * (fOld - fCur - deltaFmax) -
                        (fOld - fE) * (fOld - fE) * deltaFmax;
                    IJ.log("test = "+test);
                    if (test < 0.0) {// use new direction
                        IJ.log("new direction accepted, removing direction " +
                               maxDir);
                        for (int i=0; i < N; i++) xCur[i] = xE[i];
                        fCur = fE;
                        directions.remove(maxDir);
                        directions.add(0, newDirection);// prepend to list
                    }
                }
            }
            IJ.log("----- new round -----");
        }
        if (step >= maxNumStep)
            IJ.log("minimisation failed");
        else
            IJ.log("minimisation succeeded");
        for (int i=0; i < N; i++) fitParam[i] = xCur[i];
        return fCur;
    }

    public void showAbout() {
        //Create and set up the window.
        JFrame aboutWin = new JFrame("PlugIn Doc Browser");
        aboutWin.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        //Add scroll-pane
        JPanel p = new JPanel();
        p.setLayout(new BoxLayout(p, BoxLayout.PAGE_AXIS));
        p.setBorder(BorderFactory.createLoweredBevelBorder());
        p.setPreferredSize(new Dimension(800, 700));

        //Add small vertical space
        //p.add(Box.createRigidArea(new Dimension(0,5)));
        JScrollPane scrollPane = new JScrollPane();
        scrollPane.setAlignmentX(Component.LEFT_ALIGNMENT);
        p.add(scrollPane);

        aboutWin.add(p, BorderLayout.PAGE_START);

        // and add rest into the scrollpane

        JPanel rootpane = new JPanel();
        rootpane.setLayout(new BoxLayout(rootpane, BoxLayout.PAGE_AXIS));
        rootpane.setBorder(BorderFactory.createEmptyBorder(5, 5, 0, 5));
        
        JLabel l = new JLabel(pluginMenuName);
        Font normalFont = l.getFont();
        l.setFont(new Font(Font.SANS_SERIF, Font.BOLD,
                           normalFont.getSize()*7/5));
        l.setBorder(BorderFactory.createEmptyBorder(10, 20, 10, 20));
        rootpane.add(l);

        Font labelFont = normalFont.deriveFont(Font.ITALIC,
                                               normalFont.getSize());
        Font urlFont = new Font(Font.MONOSPACED, Font.PLAIN,
                                normalFont.getSize());

        
        rootpane.add(new JSeparator(SwingConstants.HORIZONTAL));

        java.util.Map<String,String> doc = getDocumentation();

        // extract all dictionary entries to one single page
        for (String key: doc.keySet()) {
            p = new JPanel();
            p.setLayout(new BoxLayout(p, BoxLayout.PAGE_AXIS));
            p.setBorder(BorderFactory.createEmptyBorder(5, 5, 0, 5));
            //p.setAlignmentX(Component.LEFT_ALIGNMENT);

            l = new JLabel(key);
            l.setFont(labelFont);
            //l.setAlignmentX(Component.LEFT_ALIGNMENT);
            p.add(l);
            
            JPanel p2 = new JPanel();
            p2.setLayout(new BoxLayout(p2, BoxLayout.LINE_AXIS));
            p2.add(Box.createHorizontalStrut(10));
            p2.setAlignmentX(Component.LEFT_ALIGNMENT);

            if (key.endsWith("URI") || key.endsWith("URL")) {
                try {
                    URI uri = new URI(doc.get(key));
                    JButton b = new JButton(uri.toString());
                    uris.put(b, uri);
                    b.addActionListener(this);
                    b.setFont(urlFont);
                    p2.add(b);
                }
                catch (URISyntaxException e) {
                    IJ.log("cast of \""+doc.get(key)+"\" to URI failed");
                }
            } else if (key.toLowerCase(java.util.Locale.ENGLISH).endsWith("file")) {
                URL u = getClass().getResource(doc.get(key));
                if (u != null) {
                    try {
                        URI uri = u.toURI();
                        JButton b = new JButton(uri.toString());
                        uris.put(b, uri);
                        b.addActionListener(this);
                        b.setFont(urlFont);
                        p2.add(b);
                    }
                    catch (URISyntaxException e) {
                        IJ.log("cast of \""+u+"\" to URI failed");
                    }
                } else
                    IJ.log("could not locate resource "+key+" at\n\""+
                           doc.get(key)+"\"\n");
            } else {
                JEditorPane text = new JEditorPane("text/plain", doc.get(key));
                text.setEditable(false);
                //text.setAlignmentX(Component.LEFT_ALIGNMENT);
                p2.add(text);
            }

            p.add(p2);
            p.add(Box.createVerticalStrut(5));
            p.add(new JSeparator(SwingConstants.HORIZONTAL));
            rootpane.add(p);
        }

        scrollPane.setViewportView(rootpane);

        //Display the window.
        aboutWin.pack();
        aboutWin.setVisible(true);
    }

    /** From the ActionListener interface */
    public void actionPerformed(ActionEvent e) {
        Object source = e.getSource();
        JButton button = null;
        if (source instanceof JButton) {
            button = (JButton)source;
        }
        if (button == null) return;
        if ( uris.containsKey(button) ) {
            URI uri = uris.get(source);
            try { open(uri); }
            catch (UnsupportedOperationException ex) {
                IJ.log("opening "+ uri.toString() +" failed:\n" + ex);
            }
            catch (IOException ex) {
                IJ.log("opening "+ uri.toString() +" failed:\n" + ex);
            }
        }
    }

    /** Open URI with default desktop application.
     *
     * @throws UnsupportedOperationException Either
     * java.awt.Desktop.isDesktopSupported() or desktop.isSupported(
     * java.awt.Desktop.Action.BROWSE ) returned false.
     *
     * @throws IOException The documentation of desktop.browse()
     * described under which circumstances this Exception is thrown.
     */
    public static void open(URI uri) throws UnsupportedOperationException,
                                            IOException {
        if (uri.getScheme().equals("jar")) {// need to extract from jar
            try { 
                URL url = uri.toURL();
                String suffix = 
                    (new File(uri.getSchemeSpecificPart())).getName();
                InputStream is = url.openStream();
                //IJ.log("have stream: "+is.available());
                File tmpFile = File.createTempFile("gpp",suffix);
                tmpFile.deleteOnExit();
                FileOutputStream os = new FileOutputStream(tmpFile);
                while (true) {
                    int b = is.read();
                    if (b < 0) break;
                    os.write(b);
                }
                is.close();
                uri = tmpFile.toURI();
                os.close();
            } catch (IOException e) {
                IJ.log("IOException while extracting from jar...\n"+e);
                e.printStackTrace();
            }
        }
        if (!java.awt.Desktop.isDesktopSupported())
            throw new UnsupportedOperationException("Desktop is not supported");

        java.awt.Desktop desktop = java.awt.Desktop.getDesktop();

        if (!desktop.isSupported( java.awt.Desktop.Action.BROWSE ))
            throw new UnsupportedOperationException("Desktop doesn't support the browse action");

        desktop.browse(uri);
    }


    /** Documents this plug-in using the PlugInDoc interface.
     *
     * @see PlugInDoc
     */
    public java.util.Map<String,String> getDocumentation() {
        java.util.Map<String,String> doc = 
            new java.util.LinkedHashMap<String,String>();
        doc.put("About", pluginMenuName+"is a Plugin"+
                "for liquid surface tension measurement.\n"+
    "This plug-in allows interactive adjustment of a theoretical profile\n"+
    "to an image of a pendant drop. An estimate of the quality of the fit\n"+
    "is logged to ImageJ's log window. The Plugin can also improve the fit\n"+
    "automatically by varying one or several of the parameters.");
        doc.put("Usage", 
    "Draw a rectangular ROI around the free pendant part of the drop,\n"+
    "call Plugin, check 'preview' box; then adjust parameters interactively\n"+
    "and/or fit selected parameters automatically.\n"+
    "[For more details see PDF documentation below]");
        doc.put("Author", "Adrian Daerr");
        doc.put("Version", "2013-06-06");
        doc.put("Licence", "GPL");
        doc.put("Author homepage URL",
                "http://www.msc.univ-paris-diderot.fr/~daerr/");
        doc.put("Plugin update site URL",
                "http://sites.imagej.net/Daerr/");
        doc.put("Detailed documentation PDF file", "Goutte_pendante.pdf");
        doc.put("Example image: water drop, JPEG file", "eauContrasteMax.jpg");
        return doc;
    }

}


class DropInfo {// holds information on a given ImageProcessor
    private DropInfo(ImageProcessor ip) {
        this.ip = ip;
        bounds = ip.getRoi().getBounds();
    }
    
    private final static HashMap<ImageProcessor, DropInfo> instances =
        new HashMap<ImageProcessor, DropInfo>();
    
    // Return DropInfo instance for given ImageProcessor
    public static DropInfo getDropInfo(ImageProcessor ip) {
        if (! instances.containsKey(ip)) {
            DropInfo di = new DropInfo(ip);
            instances.put(ip, di);
        }
        return instances.get(ip);
    }
    
    private ImageProcessor ip;
    public ImageProcessor getImageProcessor() { return ip; }
    
    public Rectangle bounds;
    private Area boundsArea;
    public Area getBoundsArea() {
        if (boundsArea == null)
            boundsArea = new Area(bounds);
        return boundsArea;
    }
    
    private double[] fitParam;
    public double[] getFitParam() {
        if (fitParam == null) { // idea: initialise from other ip ?
            fitParam = new double[nFitParam];
        }
        return fitParam;
    }
    public void setFitParam(double p0, double p1, double p2,
                            double p3, double p4) {
        if (fitParam == null)
            fitParam = new double[nFitParam];
        fitParam[0] = p0;
        fitParam[1] = p1;
        fitParam[2] = p2;
        fitParam[3] = p3;
        fitParam[4] = p4;
    }
    
    private double pixelScale = 1.0;
    public double getPixelScale() { return pixelScale; }
    public void setPixelScale(double s) { pixelScale = s; }
    
    private double densityGrav = 1.0;
    public double getDensityGrav() { return densityGrav; }
    public void setDensityGrav(double rho_g) { densityGrav = rho_g; }
    
    private Path2D profile;
    public Path2D getProfile() { return profile; }
    public void setProfile(Path2D profile) { this.profile = profile; }
    
    private Shape shape;
    public Shape getShape() { return shape; }
    public void setShape(Shape shape) { this.shape = shape; }
    
    public double getQuality() {
        return fitQuality(shape);
    }
        
    private double[] leftBorder;
    public double[] getLeftBorder() {
        if (leftBorder == null) findDropBorders();
        return leftBorder;
    }

    public double[] rightBorder;
    public double[] getrightBorder() {
        if (rightBorder == null) findDropBorders();
        return rightBorder;
    }

    /** Detect drop borders and store positions in the left/rightBorder arrays. */
    private void findDropBorders(ImageProcessor ip) {
        final int voisinage = 10; // pixels around border to include in fit
        final double threshold = 128.f;
        leftBorder = null;
        
        for (int y = bounds.height - 1; y >= 0; y--) {
            
            // find border positions with integer precision
            int xl = -1, xr = -1;
            for (int x = 0; x < bounds.width; x++) {
                if (ip.getPixelValue(bounds.x + x, bounds.y + y) < threshold) {
                    if (xl < 0) xl = x;  // inside drop
                    xr = -1;
                } else {// outside drop
                    if (xr < 0 && xl >= 0) xr = x;
                }
            }
            
            // don't go further if not enough pixels for subpixel-fitting 
            if (xr - xl <= voisinage ||
                xl - voisinage < 0 || xr + voisinage > bounds.width) {
                if (leftBorder != null) {
                    leftBorder[y]  = xl >= 0 ? xl - 0.5 : Double.NaN;
                    rightBorder[y] = xr >= 0 ? xr - 0.5 : Double.NaN;
                }
                continue;
            }
            
            if (leftBorder == null) { // allocate array on drop tip detection
                leftBorder = new double[y+1];
                rightBorder = new double[y+1];
            }
            
            // now determine drop borders with sub-pixel precision
            leftBorder[y]  = fitStep(xl, y, voisinage, false);
            rightBorder[y] = fitStep(xr, y, voisinage, true);
        } // end for y
    }
    
    private double fitStep(int x, int y, int n, boolean rising) {
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
}
