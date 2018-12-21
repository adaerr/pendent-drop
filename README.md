[![](https://travis-ci.com/adaerr/pendent-drop.svg?branch=master)](https://travis-ci.com/adaerr/pendent-drop)

Pendent Drop ImageJ plugin
==========================

**Questions, issues & problems**: please use the [pendent-drop topic
on the ImageJ
Forum](https://forum.image.sc/t/pendent-drop-plugin-how-to-use/290) to
ask about _installation and usage_ problems with this plugin - that
way this kind of information is easily found in one place by other
users, too, and other ImageJ users can tune in to help. The [github
issue tracker](https://github.com/adaerr/pendent-drop/issues) is meant
for _bugs and programming issues_. If in doubt go to the ImageJ Forum.

# What's the purpose of this plugin ?

_Pendent drop profile integration and fitting_

This plug-in allows for interactive or automated adjustment of a
profile to an image of a pendent drop. The surface tension, volume
and surface associated with this profile can then be obtained.

## Context

The pendent drop method is commonly used to measure surface tensions
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

For more information see the included PDF documentation (Plugins ->
Drop Analysis -> About Pendent Drop, or in the *article* directory of
the source code) or the open access publication

Daerr, A and Mogne, A 2016 *Pendent_Drop: An ImageJ Plugin to Measure
the Surface Tension from an Image of a Pendent Drop*. Journal of Open
Research Software, **4**: e3, DOI: http://dx.doi.org/10.5334/jors.97

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

# Archive

This software is archived by Zenodo:

[![DOI](https://zenodo.org/badge/18554/adaerr/pendent-drop.svg)](https://zenodo.org/badge/latestdoi/18554/adaerr/pendent-drop)
