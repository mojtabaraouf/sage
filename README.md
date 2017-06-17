<a href="http://ascl.net/1706.004"><img src="https://img.shields.io/badge/ascl-1706.004-blue.svg?colorB=262255" alt="ascl:1706.004" /></a>
[![Build Status](https://travis-ci.org/arhstevens/DarkSage.svg?branch=astevens-disc)](https://travis-ci.org/arhstevens/DarkSage)

# Dark Sage

DARK SAGE is a semi-analytic model of galaxy formation, focussed on detailing the structure and evolution of galaxies' discs.  The code-base is an extension of [SAGE](https://github.com/darrencroton/sage/) (Semi-Analytic Galaxy Evolution).  The model is described in full in the paper by [Stevens, Croton & Mutch (2016)](http://adsabs.harvard.edu/abs/2016MNRAS.461..859S).

DARK SAGE will run on any N-body simulation whose trees are organised in a supported format and contain a minimum set of basic halo properties.  Galaxy formation models built using DARK SAGE on the Millennium simulation can be downloaded at the [Theoretical Astrophysical Observatory (TAO)](https://tao.asvo.org.au/).

The code-base, written in C, should function as is, provided the required dependencies are installed.  You just need a C compiler and to point to installed GSL libraries before typing 'make'.  Once installed, please run the test script with 'python test.py' to make sure everything is working correctly.

DARK SAGE is also listed on the [Astrophysics Source Code Library](http://ascl.net/1706.004)

Queries, comments, and concerns can be emailed to Adam Stevens: astevens@swin.edu.au
