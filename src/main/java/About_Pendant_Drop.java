/*
 * Present help information about the Pendant Drop plugin
 */

import net.imagej.Dataset;
import net.imagej.ImageJ;
import org.scijava.log.LogService;
import org.scijava.ui.UIService;

import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;

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

/**
 * Show documentation and links about the Pendant Drop plugin.
 */
@Plugin(type = Command.class,
        headless = false,
        menuPath = "Plugins>Drop Analysis>About Pendant Drop")
public class About_Pendant_Drop implements Command, ActionListener {

    // needed to open images
    @Parameter
    private UIService ui;

    @Parameter
    private org.scijava.io.IOService ioService;

    @Parameter
    private LogService log;

    private final static String pluginMenuName = "Pendant drop";

    private final java.util.Map<JButton,URI> uris =
        new java.util.LinkedHashMap<JButton,URI>();

    /**
     * Show help window.
     */
    @Override
    public void run() {
        showAbout();
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

            // analyze key to decide whether we should link to some
            // resource or merely show value as plain text
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
                    log.error("cast of \""+doc.get(key)+"\" to URI failed");
                }
            } else if (key.toLowerCase(java.util.Locale.ENGLISH).endsWith("file")) {
                URL u = getClass().getResource(doc.get(key));
                if (u != null) {
                    try {
                        URI uri = u.toURI();
                        JButton b = new JButton(uri.toString());
                        uris.put(b, uri);
                        b.addActionListener(new ActionListener() {
                                public void actionPerformed(ActionEvent e) {
                                    Object source = e.getSource();
                                    if (!(source instanceof JButton)) return;
                                    JButton button = (JButton)source;
                                    if ( uris.containsKey(button) ) {
                                        URI uri = uris.get(source);
                                        try {
                                            openFile(uri);
                                        } catch (IOException x) {
                                            log.error(x);
                                        }
                                    }
                                }
                            });
                        b.setFont(urlFont);
                        p2.add(b);
                    }
                    catch (URISyntaxException e) {
                        log.error("cast of \""+u+"\" to URI failed");
                    }
                } else
                    log.error("could not locate resource "+key+" at\n\""+
                           doc.get(key)+"\"\n");
            } else if (key.toLowerCase(java.util.Locale.ENGLISH).endsWith("image")) {
                URL u = getClass().getResource(doc.get(key));
                if (u != null) {
                    try {
                        URI uri = u.toURI();
                        JButton b = new JButton(uri.toString());
                        uris.put(b, uri);
                        b.addActionListener(new ActionListener() {
                                public void actionPerformed(ActionEvent e) {
                                    Object source = e.getSource();
                                    if (!(source instanceof JButton)) return;
                                    JButton button = (JButton)source;
                                    if ( uris.containsKey(button) ) {
                                        URI uri = uris.get(source);
                                        try {
                                            openSciJava(uri);
                                        } catch (IOException x) {
                                            log.error(x);
                                        }
                                    }
                                }
                            });
                        b.setFont(urlFont);
                        p2.add(b);
                    }
                    catch (URISyntaxException e) {
                        log.error("cast of \""+u+"\" to URI failed");
                    }
                } else
                    log.error("could not locate resource "+key+" at\n\""+
                           doc.get(key)+"\"\n");

            } else {// just show as plain text
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
        if (!(source instanceof JButton)) return;
        JButton button = (JButton)source;
        if ( uris.containsKey(button) ) {
            URI uri = uris.get(source);
            try { openDesktop(uri); }
            catch (UnsupportedOperationException ex) {
                log.error("opening "+ uri.toString() +" failed:\n" + ex);
            }
            catch (IOException ex) {
                log.error("opening "+ uri.toString() +" failed:\n" + ex);
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
    public void openDesktop(URI uri) throws UnsupportedOperationException,
                                            IOException {
        if (uri.getScheme().equals("jar")) // need to extract from jar
            uri = jarURItoTmpFileURI(uri);

        if (!java.awt.Desktop.isDesktopSupported())
            throw new UnsupportedOperationException("Desktop is not supported");

        java.awt.Desktop desktop = java.awt.Desktop.getDesktop();

        if (!desktop.isSupported( java.awt.Desktop.Action.BROWSE ))
            throw new UnsupportedOperationException("Desktop doesn't support the BROWSE action");

        desktop.browse(uri);
    }

    /** Open URI using SciJava.
     *
     * @throws IOException Thrown by org.scijava.io.IOService
     */
    public void openSciJava(URI uri) throws IOException {
        if (uri.getScheme().equals("jar")) // need to extract from jar
            uri = jarURItoTmpFileURI(uri);

        // Open using appropriate plugin.
        final Dataset dataset = (Dataset)ioService.open(uri.toString());

        // display the dataset
        ui.show(dataset);
    }

    /** Open URI of a file, possibly inside a jar, with default
     * desktop application.
     *
     * @throws UnsupportedOperationException Either
     * java.awt.Desktop.isDesktopSupported() or desktop.isSupported(
     * java.awt.Desktop.Action.BROWSE ) returned false.
     *
     * @throws IOException The documentation of desktop.browse()
     * described under which circumstances this Exception is thrown.
     */
    public void openFile(URI uri) throws IOException {
        File file;
        if (uri.getScheme().equals("jar"))
            // extract from jar
            file = jarURItoTmpFile(uri);
        else
            file = new File(uri);

        if (!java.awt.Desktop.isDesktopSupported())
            throw new UnsupportedOperationException("Desktop is not supported");

        java.awt.Desktop desktop = java.awt.Desktop.getDesktop();

        if (!desktop.isSupported( java.awt.Desktop.Action.OPEN ))
            throw new UnsupportedOperationException("Desktop doesn't support the OPEN action");

        desktop.open(file);
    }

    public URI jarURItoTmpFileURI(URI uri) {
        if (!uri.getScheme().equals("jar"))
            return uri; // do nothing if not in jar
        File tmpFile = jarURItoTmpFile(uri);
        uri = tmpFile.toURI();
        return uri;
    }

    public File jarURItoTmpFile(URI uri) {
        File tmpFile = null;
        if (uri.getScheme().equals("jar")) {
            try {
                URL url = uri.toURL();
                String suffix =
                    (new File(uri.getSchemeSpecificPart())).getName();
                InputStream is = url.openStream();
                //log.debug("have stream: "+is.available());
                tmpFile = File.createTempFile("gpp",suffix);
                tmpFile.deleteOnExit();
                FileOutputStream os = new FileOutputStream(tmpFile);
                while (true) {
                    int b = is.read();
                    if (b < 0) break;
                    os.write(b);
                }
                is.close();
                os.close();
            } catch (IOException e) {
                log.error("IOException while extracting from jar...");
                log.error(e);
            }
        }
        return tmpFile;
    }

    /** Documents this plug-in using the PlugInDoc interface.
     *
     * @see PlugInDoc
     */
    public java.util.Map<String,String> getDocumentation() {
        java.util.Map<String,String> doc =
            new java.util.LinkedHashMap<String,String>();
        doc.put("About", pluginMenuName + " is a Plugin" +
                "for liquid surface tension measurement.\n" +
                "This plug-in allows for interactive or automated adjustment\n" +
                "of a theoretical profile to an image of an axisymmetric\n" +
                "pendant drop.\n" +
                "Drop properties such as surface tension, volume, interface\n" +
                "area are then calculated from the profile.");
        doc.put("Usage",
                "Draw a rectangular ROI around the free pendant part of the\n" +
                "drop before calling the Plugin; then adjust parameters\n" +
                "interactively and/or fit selected parameters automatically.\n" +
                "[For more details see PDF documentation below]");
        doc.put("Author", "Adrian Daerr");
        doc.put("Version", "2.0.0 (release date 2015-09-04)");
        doc.put("Licence", "GPL");
        doc.put("Plugin update site URL",
                "http://sites.imagej.net/Daerr/");
        doc.put("Author homepage URL",
                "http://www.msc.univ-paris-diderot.fr/~daerr/");
        doc.put("Detailed documentation PDF file",
                "article/Goutte_pendante.pdf");
        doc.put("Example image: water drop, TIFF image",
                "water_dsc1884.tif");
        return doc;
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
