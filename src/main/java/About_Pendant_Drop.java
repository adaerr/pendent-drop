/*
 * Present help information about the Pendant Drop plugin
 */

import net.imagej.ImageJ;

import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

/**
 * Show documentation and links about the Pendant Drop plugin.
 */
@Plugin(type = Command.class, headless = false, menuPath = "Plugins>Drop analysis>About Pendant Drop")
public class About_Pendant_Drop implements Command {

        @Parameter(type = ItemIO.OUTPUT)
        private String about;

        /**
         * Produce help string.
         */
        @Override
        public void run() {
                about = "Here's your documentation !";
        }

        /**
         * A {@code main()} method for testing (from IJ's Hello_World
         * example).
         * <p>
         * When developing a plugin in an Integrated Development
         * Environment (such as Eclipse or NetBeans), it is most
         * convenient to provide a simple {@code main()} method that
         * creates an ImageJ context and calls the plugin. </p>
         * <p>
         * In particular, this comes in handy when one needs to debug
         * the plugin: after setting one or more breakpoints and
         * populating the inputs (e.g. by calling something like
         * {@code ij.command().run(MyPlugin.class, "inputImage", myImage)}
         * where {@code inputImage} is the name of the field
         * specifying the input) debugging becomes a breeze.
         * </p>
         *
         * @param args unused
         */
        public static void main(final String... args) {
                // Launch ImageJ as usual.
                final ImageJ ij = net.imagej.Main.launch(args);

                // Launch our command right away.
                ij.command().run(About_Pendant_Drop.class, true);
        }

}
