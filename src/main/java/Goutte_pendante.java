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

/** An ImageJ2 plugin analyzing the shape of a pendant drop. */
@Plugin(type = Command.class,
        menuPath = "Plugins>Drop Analysis>Pendant Drop",
        initializer = "paramEstimator")
public class Goutte_pendante implements Command, Previewable {

    // -- Parameters --

    @Parameter
    private ImagePlus imp;

    @Parameter
    private LogService log;

    @Parameter(persist = false)
    private RectangleOverlay dropRegion;

    // The 'min' attribute requires a String, but
    // Double.toString(Double.MIN_VALUE) is not a constant to the
    // compiler, so we use an explicit value close to the real
    // constant
    @Parameter(persist = false, label = "tip_radius of curvature", min = "1e-300")
    private double tip_radius;

    @Parameter(persist = false, label = "capillary length", min = "1e-300")
    private double capillary_length;

    @Parameter(persist = false, label = "tip_x coordinate")
    private double tip_x;

    @Parameter(persist = false, label = "tip_y coordinate")
    private double tip_y;

    @Parameter(persist = false, label = "gravity angle (deg)")
    private double gravity_deg;

    @Parameter(persist = false, initializer = "initPixelSize")
    private double pixel_size;

    @Parameter(label = "density contrast times g")
    private double rho_g;

    @Parameter(visibility = org.scijava.ItemVisibility.MESSAGE)
    private final String label_surface_tension = "Surface tension";

    @Parameter(persist = false, visibility = org.scijava.ItemVisibility.MESSAGE, label = "Surface tension")
    private double surface_tension = 0;

    // -- Other fields --

    /** The dimensional parameters. */
    //private HashMap paramWithDim = null;

    // -- Command methods --

    @Override
    public void run() {
        log.info("drop region: +" + dropRegion.getOrigin(0)
                 + " +" + dropRegion.getOrigin(1)
                 + ", " + dropRegion.getExtent(0)
                 + " x " + dropRegion.getExtent(1));

        ImageStack stack = imp.getStack();
        // create results table

        for (int n=0; n<stack.getSize(); n++) {
            analyseImage(stack.getProcessor(n+1));
            // copy parameters into results table
        }
    }

    // -- Previewable methods --

    @Override
    public void preview() {
        updateOverlay();
        surface_tension ++;
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
    protected void paramEstimator() {
        tip_radius = 42;
        capillary_length = 43;
        tip_x = 44;
        tip_y = 45;
        gravity_deg = 0;
    }

    /** Initializes the {@link #pixel_size} parameter. */
    protected void initPixelSize() {
        pixel_size = 1;
    }

    // -- Processing --

    public void updateOverlay() {
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

        // display the dataset
        ij.ui().show(imageDisplay);

        // Launch the "CommandWithPreview" command.
        ij.command().run(Goutte_pendante.class, true);
    }

}
