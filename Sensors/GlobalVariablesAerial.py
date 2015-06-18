#-------------------------------------------------------------------------------
# Name:        Aerial SAR global variables
# Purpose:      Defines global aerial sar parameters
#                taken from Matthew Schlutz - Synthetic Aperture Radar Imaging Simulated in MATLAB
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
##PRF = 300.0 # Hz
PRF = 1200.0 # Hz
Freq = 4.5e9 # in Hz ,  Carrier frequency
Lambda = GC.C/Freq # in meters,  Carrier wavelength
H = 10000.0 # height in m
Vsat = 200.0  # in m/s
Tp = 2.5e-6 # in seconds, Range pulse length,
K = 100.0e6/Tp # Hz/s, chirp rate
La = 2.0 # in meters, antenna length

Name = "Aerial SAR Platform"
