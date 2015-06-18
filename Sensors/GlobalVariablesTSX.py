#-------------------------------------------------------------------------------
# Name:        TSX Satellite global variables
# Purpose:      Defines global TSX satellite & sensor parameters
#
# Author:      vfomin
#
# Created:     18/05/2015
# Copyright:   (c) vfomin 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# Numpy
import numpy as np

# Project
import Tools.GlobalConsts as GC

# GLOBAL VARIABLES :

Freq = 5.331e9 # in Hz ,  Carrier frequency, (http://fr.wikipedia.org/wiki/ENVISAT)
Lambda = GC.C/Freq # in meters,  Carrier wavelength
H = 785000.0 # height of the satellite in m (http://fr.wikipedia.org/wiki/ENVISAT)
Vsat = 7126.14  # in m/s (taken from here https://earth.esa.int/workshops/envisatsymposium/proceedings/posters/3P10/544447sm.pdf)
K = 4.17788e+11 # Hz/s, chirp rate (taken from here https://earth.esa.int/workshops/envisatsymposium/proceedings/posters/3P10/544447sm.pdf)
Tp = 37.120e-3 # in seconds, Range pulse length, (taken from here https://earth.esa.int/workshops/envisatsymposium/proceedings/posters/3P10/544447sm.pdf)
La = 10 # in meters, antenna length (taken from https://earth.esa.int/handbooks/asar/CNTR3-1-1.html)
PRF = 1679.9020 # Hz (taken from here https://earth.esa.int/workshops/envisatsymposium/proceedings/posters/3P10/544447sm.pdf)

Name = "Terra-Sar X"