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
        menuPath = "Plugins>OverlayNullInInit",
        initializer = "paramEstimator")
public class Goutte_pendante implements Command, Previewable {

    // -- Parameters --

    @Parameter
    private ImagePlus imp;

    @Parameter
    private LogService log;

    @Parameter(persist = false)
    private RectangleOverlay region;

    @Parameter(persist = false)
    private Overlay overlay;

    // -- Command methods --

    @Override
    public void run() {
        log.info("imp parameter in run(): "
                 + (imp==null?"":"non-") + "null");
        log.info("region parameter in run(): "
                 + (region==null?"":"non-") + "null");
        log.info("overlay parameter in run(): "
                 + (overlay==null?"":"non-") + "null");
    }

    // -- Previewable methods --

    @Override
    public void preview() {
        log.info("imp parameter in preview(): "
                 + (imp==null?"":"non-") + "null");
        log.info("region parameter in preview(): "
                 + (region==null?"":"non-") + "null");
        log.info("overlay parameter in preview(): "
                 + (overlay==null?"":"non-") + "null");
    }

    @Override
    public void cancel() {
        log.info("cancelled");
    }

    // -- Initializer methods --

    protected void paramEstimator() {
        log.info("imp parameter in initializer(): "
                 + (imp==null?"":"non-") + "null");
        log.info("region parameter in initializer(): "
                 + (region==null?"":"non-") + "null");
        log.info("overlay parameter in initializer(): "
                 + (overlay==null?"":"non-") + "null");
    }

    // -- Main method --

    /** Tests our command. */
    public static void main(final String... args) throws Exception {

        // Launch ImageJ as usual.
        final ImageJ ij = net.imagej.Main.launch(args);

        final String name = "Test Image";
        final long[] dims = { 400, 400 };
        final AxisType[] axes = { Axes.X, Axes.Y };
        final net.imagej.DatasetService datasetService = ij.dataset();
        final Dataset dataset = datasetService.create(new net.imglib2.type.numeric.integer.UnsignedByteType(), dims, name, axes);

        // create a display for the dataset
        final ImageDisplay imageDisplay =
            (ImageDisplay) ij.display().createDisplay(dataset);

        // create a rectangle
        final RectangleOverlay rectangle = new RectangleOverlay(ij.getContext());
        rectangle.setOrigin(10, 0);
        rectangle.setOrigin(20, 1);
        rectangle.setExtent(30, 0);
        rectangle.setExtent(40, 1);
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
