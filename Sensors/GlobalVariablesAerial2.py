#-------------------------------------------------------------------------------
# Name:        Aerial SAR global variables 2
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
PRF = 500.0 # Hz
Freq = 1.5e9 # in Hz ,  Carrier frequency
Lambda = GC.C/Freq # in meters,  Carrier wavelength
H = 500.0 # height in m
Vsat = 50.0  # in m/s
Tp = 2.0e-6 # in seconds, Range pulse length,
K = 1.0e14 # Hz/s, chirp rate
La = 1.5 # in meters, antenna length

Name = "Aerial SAR Platform"
