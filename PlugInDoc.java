/** <h1>PlugInDoc interface</h1>
 *
 * <p>By implementing this interface (in addition to e.g.
 * PlugIn/PlugInFilter), a plug-in can expose self-documentation in a
 * standardized way. The documentation of such plug-in can then be
 * retrieved or displayed by other tools (e.g. using the
 * Show_help_for_PlugIn plugin from within ImageJ).</p>
 *
 * <p>The aim is to make it easy for plug-in-writers to point to whatever
 * documentation they wish to make available (brief text,
 * javadoc-generated API, ...), thus making it easier for users to
 * find this documentation.</p>
 *
 * <p>The single method returns a dictionary of individual 'help
 * items'. Any String is allowed as a key or value, we however suggest
 * that the key be a human readable description of the meaning of the
 * value (e.g. "Author", "Licence", "Algorithm description"). In
 * addition, we suggest that a key should end in "URI" or "URL" to
 * indicate that the corresponding value may be a Universal Resource
 * Identifier. This allows to link freely to ressources on the WWW, in
 * jar-files etc.</p>
 *
 * <h2>Example usage 1</h2>
 *  <pre>
 * import ij.plugin.PlugIn;
 * import java.util.Map;
 * [...]
 *
 * public class Get_help_for_PlugIn implements PlugIn, PlugInDoc {
 *
 *     public void run(String arg) {
 *         [...]                     // the usual PlugIn stuff here
 *     }
 *
 *     // and here is the single method required by PlugInDoc:
 *
 *     public Map&lt;String,String> getDocumentation() {
 *         // using a LinkedHashMap is recommended to preserve enumeration order
 *         Map&lt;String,String> doc = new LinkedHashMap&lt;String,String>();
 *
 *         // simple, compact documentation (compare to other example):
 *
 *         doc.put("About", "Get_help_for_PlugIn\n"+
 *                          "\n"+
 *                          "Extract and show documentation about plug-ins.\n"+
 *                          " If the description field ends in \"URI\" or \"URL\", \n"+
 *                          "and the value can indeed be parsed as a URI, offers \n"+
 *                          "to open the link with the system default application");
 *                          "\n"+
 *                          "Adrian Daerr\n"+
 *                          "Paris 2009\n"+
 *                          "public domain");
 *
 *         doc.put("Homepage URL",
 *                          "http://www.msc.univ-paris-diderot.fr/~daerr/");
 *
 *         return doc;
 *     }
 * }
 * </pre>
 *
 * <h2>Example usage 2  (see also the <a href="http://www.msc.univ-paris-diderot.fr/~daerr/ijstuff/source/Get_doc_for_PlugIn.html">source of Get_doc_for_PlugIn</a>)</h2>
 *
 * <pre>
 *     // [rest of plug-in not shown]
 *     public Map&lt;String,String> getDocumentation() {
 *         Map&lt;String,String> doc = new LinkedHashMap&lt;String,String>();
 *
 *         doc.put("About",   "Get_help_for_PlugIn\n\n"+
 *                            "Extract and show documentation about plug-ins.\n"+
 *                            " If the description field ends in \"URI\" or \"URL\", \n"+
 *                            "and the value can indeed be parsed as a URI, offers \n"+
 *                            "to open the link with the system default application");
 *         doc.put("Author",  "Adrian Daerr");
 *         doc.put("Version", "2009-04-14");
 *         doc.put("Licence", "public domain");
 *         doc.put("API URL", "http://www.msc.univ-paris-diderot.fr/~daerr/ijstuff/api/Get_help_for_PlugIn.html");
 *         doc.put("Browsable source code URL",
 *                            "http://www.msc.univ-paris-diderot.fr/~daerr/ijstuff/source/Get_help_for_PlugIn.html");
 *         doc.put("Class update URL",
 *                            "http://www.msc.univ-paris-diderot.fr/~daerr/ijstuff/source/Get_help_for_PlugIn.class");
 *         doc.put("Source update URL",
 *                            "http://www.msc.univ-paris-diderot.fr/~daerr/ijstuff/source/Get_help_for_PlugIn.java");
 *         doc.put("Author homepage URL",
 *                            "http://www.msc.univ-paris-diderot.fr/~daerr/");
 *         doc.put("E-mail URI", "mailto:author@does-not-exist.imaginary");
 *         return doc;
 *     }
 * </pre>
 */

public interface PlugInDoc {
    /** Documentation Strings in a Map. See introduction to this
     * interface for a more detailed description and implementation
     * examples. */
    public java.util.Map<String,String> getDocumentation();
}
