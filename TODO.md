List of improvements for the Pendent Drop plugin
-------------------------------------------------------------------

- calculate and output information that helps assessing the reliability of the measurement
 * Bond number (ratio of weight to maximum capillary retaining force), indicates how far we are from break-off (thus how much the shape is influenced by gravity)
 * error bars from Hessian matrix around optimum parameter set
- more precise/informative display of parameters
 * currently 0.0012345 is shown to 3 decimal places only, i.e. 0.001
 * initial parameters are shown with an unreasonable amount of (in)significant digits
 * parameter change buttons modify by units, which rarely makes sense (+-1% or so instead ?)
 * change in scale changes drop contour: at least when previous value was 1 (arguably equivalent to 'no scale'), we should probably rescale the other variables to keep the visual contour in place.
- rename 'pixel size' parameter for more clarity
- add possibility to output profile slope at specific points (e.g. on selection boundary), for contact angle measurement
- add possibility to fit a sessile drop and reflection simultaneously (still with optional contact angle output)

- catch case where plugin is called without open image
  ==? check that display!=null in paramInitializer
  (cf https://forum.image.sc/t/pendent-drop-plugin-how-to-use/290/31 )
- reorganise/clean source code
 * split into separate classes
 * separate optimisation code
- provide alternative optimisation algorithms (e.g. quadratic descent
  using Hessian, mixed strategies)
- provide possibility of parallel calculation on stacks
