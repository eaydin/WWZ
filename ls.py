#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ls.py
# 16 May 22:31
# 11 September 12:16 semicolon fixes
# Copyright 2014 M. Emre Aydin <emre.m.aydin@gmail.com>

"""
    A Simple Lomb-Scargle (LS) Code to Check the Research for the
    Weighted Wavelet Z-Transform Translation of Python.
    The script uses the LS method in the SciPy Package.
    Available at http://github.com/eaydin
"""

# Checking the necessary modules and importing them if exists
import sys
try:
    import numpy as np
except:
    print "Please Install Numpy"
    raise SystemExit

try:
    import scipy as sp
except:
    print "Please Install SciPy"
    raise SystemExit

try:
    import argparse
except:
    "Argparse Not Available?! You sure this is Python 2.7?"
    raise SystemExit

from scipy.signal import spectral

class LombScargle(object):
    """The Main class object.
    This is the main class object to create a Lomb-Scargle periodogram.
    Usage:
    LombScargle(timeData, magnitudeData, freq_low, freq_high, \
                freq_num)
    
    Arguments:
    timeData: float array of time values in JD
    magnitudeData: float array of magnitude values
    freq_low: the low frequency value
    freq_high: the high frequency value
    freq_num: number of total frequencies to evaluate
    
    Available Attributes:
    'power': returns the LS-Power only.
    'freqs': returns the frequencies calculated.
    'periodogram': returns the LS Periodogram as [freqs, power].
    '    
    """
    
    def __init__(self, timeData, magnitudeData, freq_low, \
                 freq_high, freq_num):
        """Initializing the LombScargle Object"""               
        self.freqs = np.linspace(freq_low, freq_high, freq_num)
             
        # Convert arrays to Numpy Arrays for lombscargle to understand
        timeArray = np.asarray(timeData)
        magnitudeArray = np.asarray(magnitudeData)
        
        # Calculating the periodogram and return it
        self.power = spectral.lombscargle(timeArray, magnitudeArray, \
                                          self.freqs)
        self.periodogram = np.column_stack([self.freqs, self.power])
        
        # The Frequency and LS-Power of the Maximum LS-Power
        self.maxLS = [self.freqs[self.power.argmax(axis=0)], \
                      np.amax(self.power)]
        
        
if __name__ == '__main__':
    """The main function run when standalone"""
    
    description = """
    The Lomb-Scargle Periodogram Program.
    This script is written for research purposes to accompany the 
    Weighted Wavelet Z-Transform implementation.
    M. Emre Aydin - emre.m.aydin@gmail.com
    http://about.me/emre.aydin
    Available at http://github.com/eaydin
    
    Import this script via Python to use it as a module, rather than
    a standalone script. (import ls)    
    """
    
    # Parsing arguments
    
    parser = argparse.ArgumentParser(prog="ls.py", \
                    formatter_class=argparse.RawDescriptionHelpFormatter, \
                    description=description)
    parser.add_argument("-f", "--file", type=argparse.FileType("r"), \
                        default=sys.stdin, required=True, \
                        help="the Input file, Raw Lightcurve")
    parser.add_argument("-o", "--output", type=argparse.FileType("w"), \
                        help="the Output file name")
    parser.add_argument("-s", "--stdout", action="store_true", \
                        default=False, help="Prints the results to \
                        the standard output. Automatically set to 'True' \
                        if '-o' option is not defined, otherwise 'False'.")
    parser.add_argument("-l", "--freq-low", type=float, required=True, \
                        help="the Low Frequency Value")
    parser.add_argument("-hi", "--freq-high", type=float, required=True, \
                        help="the High Frequency Value")
    parser.add_argument("-n", "--freq-num", type=int, required=False, \
                        default=1000, help="the Number of Frequencies to \
                        calculate. Input type is integer. Default value \
                        is 1000.")
    parser.add_argument("--no-headers", action="store_true", \
                        default=False, help="Doesn't print headers to output \
                        if set. Default is 'False'.")
    args = parser.parse_args()

    # Check if not output file is specified, then set stdout to 'True'
    if args.output is None:
        args.stdout = True
    
    magnitudeData = []
    timeData = []

    # Read the lightcurve data
    # Dismiss lines beginning with '%' and '#'    
    for line in args.file:
        if line.strip()[0] != "%" and line.strip()[0] != "#":
            timeData.append(float(line.split()[0].strip()))
            magnitudeData.append(float(line.split()[1].strip()))
      
    s = LombScargle(timeData, magnitudeData, args.freq_low, \
                    args.freq_high, args.freq_num)   

    # Setting the float precision of Numpy to 5
    np.set_printoptions(precision=5)
    # Suppress the printing options for Numpy
    np.set_printoptions(suppress=True)
    # Print the full Numpy Array
    np.set_printoptions(threshold='nan')

    if args.stdout:
        if not args.no_headers:
            print "# Frequency \t LS-Power"
        for i in s.periodogram:
            print "%.5f %.5f "% (i[0], i[1])

    if args.output is not None:
        if not args.no_headers:
            np.savetxt(args.output, s.periodogram, delimiter="\t", \
                       fmt="%10.5f")
        else:
            np.savetxt(args.output, s.periodogram, fmt="%10.5f", \
                       delimiter="\t", comments="#", \
                       header=" Frequency \t LS-Power")
