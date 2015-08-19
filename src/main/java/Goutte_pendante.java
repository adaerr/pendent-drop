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

import ij.ImagePlus;
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
        menuPath = "Plugins>Drop Analysis>Pendant Drop")
public class Goutte_pendante implements Command, Previewable {

    // -- Parameters --

    @Parameter
    private ImagePlus imp;

    @Parameter(persist = false, initializer = "initTitle")
    private String title;

    @Parameter
    private LogService log;

    // -- Other fields --

    /** The original title of the image. */
    private String initialTitle;

    // -- Command methods --

    @Override
    public void run() {
        // Set the image's title to the specified value.
        imp.setTitle(title);
        log.info("roi: "+imp.getProcessor().getRoi());
    }

    // -- Previewable methods --

    @Override
    public void preview() {
        run();
    }

    @Override
    public void cancel() {
        // Set the image's title back to the original value.
        imp.setTitle(initialTitle);
    }

    // -- Initializer methods --

    /** Initializes the {@link #title} parameter. */
    protected void initTitle() {
        title = initialTitle = imp.getTitle();
    }

    // -- Main method --

    /** Tests our command. */
    public static void main(final String... args) throws Exception {
        final String testImagePath = "/home/adrian/Programmes/plugins_ImageJ_src/Traitement_Gouttes/src/main/resources/article/eauContrasteMax.jpg";

        // Launch ImageJ as usual.
        //final ImageJ ij = net.imagej.Main.launch(args);
        final ImageJ ij = new ImageJ();
        //ij.ui().showUI();
        ij.ui().showUI("swing");

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
