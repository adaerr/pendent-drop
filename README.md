Pendant Drop ImageJ plugin
==========================

# Pendant drop profile integration and fitting

This plug-in allows for interactive or automated adjustment of a
profile to an image of a pendant drop. The surface tension, volume
and surface associated with this profile can then be obtained.

## Context

The pendant drop method is commonly used to measure surface tensions
of liquids. It consists in analysing the shape of a drop hanging
typically from a capillary tube and about to detach (sometimes the
inverse situation of a bubble forming at the bottom of a liquid is
preferred, or that of a sessile drop or bubble). The shape is very
sensitive to the unknown interfacial tension. The drop profile is
described by only one non-dimensional parameter (tip radius over
capillary length), although in practice five dimensional parameters
can be adjusted within this plug-in: tip position and curvature, tilt
of symetry axis and capillary length. The surface tension is
calculated from the latter if the density difference is given.

For more information see the included PDF documentation
(Plugins -> Drop Analysis -> About Pendant Drop)

# Installation

The latest stable version is uploaded to the Fiji update site
http://sites.imagej.net/Daerr/

Follow the procedure described on the Fiji wiki at
http://fiji.sc/How_to_follow_a_3rd_party_update_site
to add that site to Fiji's list and install the plugin.

I do not recommend using the git version, simply because it may be in
an unusable state for long periods during rewrites. If you
nevertheless insist on trying the latest version under development,
check the commit history and use either maven or fiji+jar to compile
and package the plugin for you.

# Author, copyright, distribution policy and disclaimer

Written by Adrian Daerr

Copyright Adrian Daerr & Université Paris Diderot

This plugin is distributed under the GNU General Public Licence,
version 3. The licencing text can be found on the site of the Free
Software Foundation at: http://fsf.org/

If you require another licence please contact the author or
his employer.

THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
APPLICABLE LAW. THE PROGRAM IS PROVIDED "AS IS" WITHOUT WARRANTY OF
ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE
PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME
THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
