#!/usr/bin/env python


"""
    A Weighted Wavelet Z-Transformation Application for Python.
    Translated by M. Emre Aydin - emre.m.aydin@gmail.com
    http://about.me/emre.aydin
    Available at https://github.com/eaydin/WWZ

For details, read the docstrings of WWZ class.
"""

import math
import sys
import os
from datetime import datetime

# Checking if Python is 2.7
if not sys.version_info[:2] == (2, 7):
    print("You need Python 2.7 for this script to run.")
    raise SystemExit

# Checking if we have numpy
try:
    import numpy as np
except:
    print("Please Install Numpy")
    raise SystemExit

# Checking if we have argparse (which is default in 2.7 so this is
# probably a useless check since we already check if Python is 2.7)
try:
    import argparse
except:
    print("Argparse Not Available?! You sure this is Python 2.7?")
    raise SystemExit


# The WWZ Class Begins

class WWZ(object):
    """The Main class object.
    This object class does not get any arguments.
    Available methods are:
        readfile()
        roundtau()
        maketau()
        makefreq()
        matrix_inv()
        wwt()
        writefile()
        writegnu()

    Arguments:
        fileName: the input filename, should be the lightcurve.
        outputfileName: the output filename.
        flo: the low frequency value. Float.
        fhi: the high frequency value. Float.
        df: the frequency step. Float.
        dcon: the C Window constant. Float.
        timedivisions: the The Divisions value, choose 50.0 if not
            sure, that is the default used by Templeton.
        max_periods: set True if you want the second output file.
            It outputs the tau's with maximum period values for easier
            estimation.
        gnuplot_compatible: splits tau values by a blank line if set
            True, so that pm3d of gnuplot easily maps the plot.
    """

    def __init__(self):
        """Initializing the object"""

    def readfile(self, fileName):
        """Read the input file.
        The argument is the file pointer, not the filename as a string.
        The values in file should be delimited with spaces or tabs.
        Ignores lines starting with # and %, as if they're comment lines.

        Returns two arrays:
            Time value, read from the first column of input file.
            Magnitude value, read from the second column of input file.
        """

        time = []
        magnitude = []

        for line in fileName:
            # Check if it's a comment line
            if line.strip()[0] != "%" and line.strip()[0] != "#":
                line_time = float(line.split()[0])
                line_mag = float(line.split()[1])
                time.append(line_time)
                magnitude.append(line_mag)

        fileName.close()

        # Just a routine check for parameter number equality
        # This should be cleaned up a bit
        if len(time) != len(magnitude):
            print("Number of Time and Magnitude input do not match. \
                  Please check the input file.")
            raise SystemExit

        # Return two arrays
        return time, magnitude

    def roundtau(self, darg):
        """Rounds the tau's. from G. Foster's Code.
        This is actually called by the maketau method.
        The input is dtspan/timedivisions,
        where dtspan is the entire timespan of the lightcurve.
        so dtspan = time[-1] - time[0]

        Returns the round value.
        """

        dex = math.pow(10, int(math.log(darg, 10)))

        darg = darg / dex

        if darg >= 5:
            darg = 5.0
        elif darg >= 2:
            darg = 2.0
        else:
            darg = 1.0

        darg = darg * dex
        return darg


    def maketau(self, time, timedivisions):
        """The maketau method.
        Arguments are:
            time = The array of time values
            timedivisions = The value of timedivisions to create tau values

            Returns an array of calculated tau values.
        """
        # The Maketau section
        # Lines 90 - 122 Fortran

        dtauhi = time[-1]
        dtaulo = time[0]

        dtspan = dtauhi - dtaulo
        dtstep = roundtau(dtspan / timedivisions)

        dtaulo = dtstep * int(dtaulo / dtstep)
        dtauhi = dtstep * int((dtauhi / dtstep) + 0.5)

        return np.arange(dtaulo, dtauhi + dtstep, dtstep).tolist()

    def makefreq(self, flo, fhi, df):
        """The Makefreq section.

        Arguments are:
            flo = Low Frequency
            fhi = High Frequency
            df = Frequency Step
        """
        # Lines 149 - 181 Fortran
        # Lines 356 - 370 Java
        nfreq = int((fhi - flo) / df) + 1
        return [(flo + ((i - 1) * df)) for i in range(1, nfreq + 1)]

    def wwt(self, time, magnitude, flo, fhi, df, dcon, timedivisions):
        """The WWZ Algorithm
        Arguments are:
            time = The time values as an array
            magnitude = The magnitude values as an array
            flo = The Low Frequency
            fhi = The High Frequency
            df = The Frequency Step
            dcon = The C constant of WWZ Window
            timedivisions = The TAU steps

            Returns a NumPy Array
        """

        dave = np.mean(magnitude)
        dvar = np.var(magnitude)
        
        freq = self.makefreq(flo, fhi, df)
        nfreq = len(freq)
        dmat = np.zeros(shape=(3,3))

        ### End of Initializing

        tau = self.maketau(time, timedivisions)
        ntau = len(tau)

        ### WWT Stars Here

        dvec = [0,0,0] # length is 3
        dcoef = [0,0,0] # length is 3

        itau = 0
        ifreq = 0
        idat = 0

        domega = 0.0
        dweight2 = 0.0
        dz = 0.0
        dweight = 0.0

        dcc = 0.0
        dcw = 0.0
        dss = 0.0
        dsw = 0.0
        dxw = 0.0
        dvarw = 0.0

        dtau = 0.0

        dpower = 0.0
        dpowz = 0.0
        damp = 0.0
        dneff = 0.0
        davew = 0.0

        dfre = 0.0

        n1 = 0
        n2 = 0

        dmz = 0.0
        dmzfre = 0.0
        dmzamp = 0.0
        dmcon = 0.0
        dmneff = 0.0

        twopi = 2.0 * math.pi

        ndim = 2
        itau1 = 0 # --> 1 or 0 ??
        itau2 = ntau # --> ???
        ifreq1 = 1
        ifreq2 = nfreq
        nstart = 1

        # Creating output arrays
        output = np.empty((ntau*(nfreq-1), 6))
        numdat=len(time)
        index = 0

        # Use for itau in range(itau1,itau2) for parallel
        for itau in range(0, itau2):
            nstart = 1
            dtau = tau[itau]

            dmfre = 0.0
            dmamp = 0.0
            dmcon = 0.0
            dmneff = 0.0
            dmz = -1.0 # less than the smallest WWZ
            
            for ifreq in range(ifreq1, ifreq2):
                dfre = freq[ifreq]
                domega = dfre * twopi

                for i in range(0, ndim + 1):
                    dvec[i] = 0.0
                    for j in range(0, ndim + 1):
                        dmat[i][j] = 0.0

                dweight2 = 0.0

                for idat in range(nstart, numdat):

                    dz = domega * (time[idat] - dtau)
                    dweight = math.exp(-1.0 * dcon * dz * dz)

                    if (dweight > 10**(-9)):
                        dcc = math.cos(dz)
                        dcw = dweight * dcc
                        dss = math.sin(dz)
                        dsw = dweight * dss
                        dmat[0][0] = dmat[0][0] + dweight
                        dweight2 = dweight2 + (dweight**2)
                        dmat[0][1] = dmat[0][1] + dcw
                        dmat[0][2] = dmat[0][2] + dsw
                        dmat[1][1] = dmat[1][1] + (dcw * dcc)
                        dmat[1][2] = dmat[1][2] + (dcw * dss)
                        dmat[2][2] = dmat[2][2] + (dsw * dss)

                        dxw = dweight * magnitude[idat]
                        dvec[0] = dvec[0] + dxw
                        dvarw = dvarw + (dxw * magnitude[idat])
                        dvec[1] = dvec[1] + (dcw * magnitude[idat])
                        dvec[2] = dvec[2] + (dsw * magnitude[idat])

                    elif dz > 0.0:
                        break
                    else:
                        nstart = idat + 1

                dpower = 0.0
                damp = 0.0

                for n1 in range(0, ndim + 1):
                    dcoef[n1] = 0.0

                if (dweight2 > 0.0):
                    dneff = (dmat[0][0] * dmat[0][0]) / dweight2
                else:
                    dneff = 0.0

                if (dneff > 3.0):

                    for n1 in range(0, ndim + 1):
                        dvec[n1] = dvec[n1] / dmat[0][0]

                        for n2 in range(1, ndim + 1 ):
                            dmat[n1][n2] = dmat[n1][n2] / \
                                                dmat[0][0]
                
                    if (dmat[0][0] > 0.005):
                        dvarw = dvarw / dmat[0][0]
                    else:
                        dvarw = 0.0
                    

                    dmat[0][0] = 1.0
                    davew = dvec[0]
                    dvarw = dvarw - (davew ** 2)

                    if (dvarw <= 0.0):
                        dvarw = 10**-12

                    for n1 in range(1, ndim + 1):
                        for n2 in range(0, n1):
                            dmat[n1][n2] = dmat[n2][n1]

                    dmat = np.linalg.inv(dmat)

                    for n1 in range(0, ndim + 1):
                        for n2 in range(0, ndim + 1):
                            dcoef[n1] = dcoef[n1] + dmat[n1][n2] * \
                                                    dvec[n2]

                        dpower = dpower + (dcoef[n1] * dvec[n1])

                    dpower = dpower - (davew ** 2)
                    dpowz = (dneff - 3.0) * dpower / (dvarw - dpower) / 2.0
                    dpower = (dneff - 1.0) * dpower / dvarw / 2.0
                    damp = math.sqrt(dcoef[1] * dcoef[1] + \
                                     dcoef[2] * dcoef[2])

                else:

                    dpowz = 0.0
                    dpower = 0.0
                    damp = 0.0

                    if (dneff < (10**(-9))):
                        dneff = 0.0

                if (damp < (10**(-9))):
                    damp = 0.0

                if (dpower < (10**(-9))):
                    dpower = 0.0

                if (dpowz < (10**(-9))):
                    dpowz = 0.0

                # Let's write everything out.
                output[index] = [dtau, dfre, dpowz, damp, dcoef[0], dneff]

                index = index + 1

                if (dpowz > dmz):
                    dmz = dpowz
                    dmfre = dfre
                    dmamp = damp
                    dmcon = dcoef[0]
                    dmneff = dneff

        return output  


    def writefile(self, wwz_output, outputFile, no_headers, max_periods):
        """The write file method.
        Arguments are:
            wwz_output = The NumPy array of WWZ values to write
            outputFile = The output file pointer, not the filename
            no_headers = If true, will not write headers to the output
            max_periods = If true, will create a file with period values
                          of maximum WWZ statistics
        """

        np.set_printoptions(precision=5)
        np.set_printoptions(suppress=True)
        np.set_printoptions(threshold='nan')

        if no_headers:
            np.savetxt(outputFile, wwz_output, delimiter="\t", \
                          fmt="%10.4f")
        else:
            np.savetxt(outputFile, wwz_output, delimiter="\t", \
                        fmt="%10.4f", comments="#", \
                        header="%9s %10s %10s %10s %10s %10s" % \
                        ("TAU", "FREQ", "WWZ", "AMP", "COEF", "NEFF"))

    def writegnu(self, wwz_output, outputFile, no_headers, \
                max_periods, ntau):
        """The write file method, adapted to work with GnuPlot.
        Arguments are:
            wwz_output = The NumPy array of WWZ values to write
            outputFile = The output file pointer, not the filename
            no_headers = If true, will not write headers to the output
            max_periods = If true, will create a file with period values
                          of maximum WWZ statistics
            ntau = The number of tau values, this is needed in order to
                   split the wwz_output equally
                To calculate ntau, use the equation below:
                    len(wwz_output) /
                      int(((freq_high - freq_low) / freq_step) + 1)
        """

        np.set_printoptions(precision=5)
        np.set_printoptions(suppress=True)
        np.set_printoptions(threshold='nan')

        splitArray = np.vsplit(wwz_output, ntau)

        # check if the file is in append mode
        # if not, reopen it
        if outputFile.mode != 'a':
            outputFile.close()
            outputFile = open(outputFile.name, "a")

        # write headers if expected
        if not no_headers:
            header="%9s %10s %10s %10s %10s %10s" % \
                   ("TAU", "FREQ", "WWZ", "AMP", "COEF", "NEFF")
            outputFile.write("#" + header + "\n")

        # split the array and add newlines in between tau values,
        # write the output
        for i in range(0,ntau):
            np.savetxt(outputFile, splitArray[i], delimiter="\t", \
                          fmt="%10.4f")
            if i != ntau-1:
                outputFile.write("\n")

# The WWZPAR Class Begins

class WWZPAR(object):
    """This is for Parallel Processing only.
    This works different than the WWZ class, it takes arguments directly.

    It gets input as a filepointer, NOT as arrays!
    Arguments:
        fileName: the input filename, should be the lightcurve.
        outputfileName: the output filename.
        flo: the low frequency value. Float.
        fhi: the high frequency value. Float.
        df: the frequency step. Float.
        dcon: the C Window constant. Float.
        timedivisions: the The Divisions value, choose 50.0 if not sure,
                        that is the default used by Templeton.
        max_periods: set True if you want the second output file.
                      It outputs the tau's with maximum period values
                      for easier estimation.
        gnuplot_compatible: splits tau values by a blank line if set True,
                             so that pm3d of gnuplot easily maps the plot.
    """

    def __init__(self, fileName, outputfileName, flo, fhi, df, dcon, \
                timedivisions, max_periods, gnuplot_compatible):
        """Initializing the object"""

        self.inputfile = fileName

        self.outputfilename1 = outputfileName
        self.outputfilename2 = outputfileName.name + ".max_periods"

        self.max_periods = max_periods
        self.gnuplot_compatible = gnuplot_compatible

        self.timedivisions = timedivisions
        # This (50) is an assumption by Templeton.
        # VStars leaves this optional
        # but keeps the default value

        self.fhi = fhi
        self.flo = flo
        self.df = df

        # dcon is the Window Constant "c" in Foster's equations.
        self.dcon = dcon

        self.fileName = fileName

        self.time = [] # Input time, first column in the file
        self.magnitude = [] # Input magnitude, second column in the file

        self.dave = 0.0 # average
        self.dvar = 0.0 # variance

        self.nfreq = int((self.fhi - self.flo) / self.df) + 1

        self.freq = []

        self.dmat = np.zeros(shape=(3, 3))


    def readfile(self):
        """Read the input file"""

        #read_file = open(self.inputfilename, "r")
        for line in self.inputfile:
            if line.strip()[0] != "%" and line.strip()[0] != "#":
                line_time = float(line.split()[0])
                line_mag = float(line.split()[1])
                self.time.append(line_time)
                self.magnitude.append(line_mag)
                self.dave = self.dave + line_mag
                self.dvar = self.dvar + (line_mag ** 2)

        self.inputfile.close()
        if len(self.time) != len(self.magnitude):
            print("Number of Time and Magnitude input do not match. \
                   Please check the input file.")
            raise SystemExit
            # Just a routine check for parameter number equality.
            # This should be cleaned up a bit.

        # Calculating Header Values
        self.numdat = len(self.time)
        self.dave = self.dave / self.numdat
        self.dvar = (self.dvar / self.numdat) - (self.dave ** 2)
        self.dsig = math.sqrt((self.dvar * self.numdat) / (self.numdat - 1))


    def roundtau(self, darg):
        """Rounds the tau's. from G. Foster's Code."""

        dex = math.log(darg, 10)
        nex = int(dex)

        darg = darg / math.pow(10, nex)

        if darg >= 5:
            darg = 5.0
        elif darg >= 2:
            darg = 2.0
        else:
            darg = 1.0

        darg = darg * math.pow(10, nex)
        return darg

    def matrix_inv(self,input_matrix):
        """The Matrix Inversion Function"""
        # Lines 202 - 252 Fortran

        ndim = 2
        dsol = np.zeros(shape=(3, 3))

        for i in range(0, 3):
            for j in range(0, 3):
                dsol[i][j] = 0.0
            dsol[i][i] = 1.0

        for i in range(0, 3):
            if input_matrix[i][i] == 0.0:
                if i == ndim:
                    return
                for j in range(0, 3):
                    if input_matrix[j][i] != 0.0:
                        for k in range(0, 3):
                            input_matrix[i][k] = input_matrix[i][k] + \
                                                 input_matrix[j][k]
                            dsol[i][j] = dsol[i][j] + dsol[j][k]

            dfac = input_matrix[i][i]

            for j in range(0, 3):
                input_matrix[i][j] = input_matrix[i][j] / dfac
                dsol[i][j] = dsol[i][j] / dfac

            for j in range(0, 3):
                if j != i:
                    dfac = input_matrix[j][i]
                    for k in range(0, 3):
                        input_matrix[j][k] = input_matrix[j][k] - \
                                            (input_matrix[i][k] * dfac)
                        dsol[j][k] = dsol[j][k] - (dsol[i][k] * dfac)

        return dsol


    def maketau(self):
        """The Maketau section"""
        # Lines 90 - 122 Fortran

        dtaulo = self.time[0] # the java translation uses 1 for this, 
        # but that should be a bug.
        dtauhi = self.time[-1]

        dtspan = dtauhi - dtaulo
        dtstep = self.roundtau(dtspan / self.timedivisions)

        dtaulo = dtstep * int(dtaulo / dtstep)
        dtauhi = dtstep * int((dtauhi / dtstep) + 0.5)

        self.tau = []

        dtau = dtaulo
        while dtau <= dtauhi:
            self.tau.append(dtau)
            dtau = dtau + dtstep


        self.ntau = len(self.tau)

    def makefreq(self):
        """The Makefreq section"""
        # Lines 149 - 181 Fortran
        # Lines 356 - 370 Java

        self.freq.append(self.flo)

        # These lines seem skeptical!
        for i in range(1,self.nfreq+1):
            self.freq.append(self.flo + ((i - 1) * self.df))

    def wwt(self, output1_par, itau1_par, itau2_par):
        """The WWZ Algorithm in Parallel Mode"""
        
        output1 = open(output1_par, "w")
        max_periods = self.max_periods
        gnuplot_compatible = self.gnuplot_compatible

        if max_periods:
            output2_par = output1_par + ".max_periods.par"
            output2 = open(output2_par, "w")

        dvec = [0,0,0] # length is 3
        dcoef = [0,0,0] # length is 3

        itau = 0
        ifreq = 0
        idat = 0

        domega = 0.0
        dweight2 = 0.0
        dz = 0.0
        dweight = 0.0

        dcc = 0.0
        dcw = 0.0
        dss = 0.0
        dsw = 0.0
        dxw = 0.0
        dvarw = 0.0

        dtau = 0.0

        dpower = 0.0
        dpowz = 0.0
        damp = 0.0
        dneff = 0.0
        davew = 0.0

        dfre = 0.0

        n1 = 0
        n2 = 0

        dmz = 0.0
        dmzfre = 0.0
        dmzamp = 0.0
        dmcon = 0.0
        dmneff = 0.0

        twopi = 2.0 * math.pi

        ndim = 2
        itau1 = 0 # ----> 1 or 0 ??
        itau2 = self.ntau # -------> ???
        ifreq1 = 1
        ifreq2 = self.nfreq
        nstart = 1

        dmat_par = np.zeros(shape=(3, 3))

        for itau in range(itau1_par, itau2_par):
            nstart = 1
            dtau = self.tau[itau]

            dmfre = 0.0
            dmamp = 0.0
            dmcon = 0.0
            dmneff = 0.0
            dmz = -1.0 # less than the smallest WWZ

            for ifreq in range(ifreq1, ifreq2 + 1):
                dfre = self.freq[ifreq]
                domega = dfre * twopi

                for i in range(0, ndim + 1):
                    dvec[i] = 0.0
                    for j in range(0, ndim + 1):
                        dmat_par[i][j] = 0.0

                dweight2 = 0.0

                for idat in range(nstart, self.numdat):

                    dz = domega * (self.time[idat] - dtau)
                    dweight = math.exp(-1.0 * self.dcon * dz * dz)

                    if dweight > 10**(-9):
                        dcc = math.cos(dz)
                        dcw = dweight * dcc
                        dss = math.sin(dz)
                        dsw = dweight * dss
                        dmat_par[0][0] = dmat_par[0][0] + dweight
                        dweight2 = dweight2 + (dweight**2)
                        dmat_par[0][1] = dmat_par[0][1] + dcw
                        dmat_par[0][2] = dmat_par[0][2] + dsw
                        dmat_par[1][1] = dmat_par[1][1] + (dcw * dcc)
                        dmat_par[1][2] = dmat_par[1][2] + (dcw * dss)
                        dmat_par[2][2] = dmat_par[2][2] + (dsw * dss)

                        dxw = dweight * self.magnitude[idat]
                        dvec[0] = dvec[0] + dxw
                        dvarw = dvarw + (dxw * self.magnitude[idat])
                        dvec[1] = dvec[1] + (dcw * self.magnitude[idat])
                        dvec[2] = dvec[2] + (dsw * self.magnitude[idat])

                    elif dz > 0.0:
                        break
                    else:
                        nstart = idat + 1

                dpower = 0.0
                damp = 0.0

                for n1 in range(0, ndim + 1):
                    dcoef[n1] = 0.0

                if dweight2 > 0.0:
                    dneff = (dmat_par[0][0] * dmat_par[0][0]) / dweight2
                else:
                    dneff = 0.0

                if dneff > 3.0:

                    for n1 in range(0, ndim + 1):
                        dvec[n1] = dvec[n1] / dmat_par[0][0]

                        for n2 in range(1, ndim + 1 ):
                            dmat_par[n1][n2] = dmat_par[n1][n2] / \
                                               dmat_par[0][0]

                    if dmat_par[0][0] > 0.0:
                        dvarw = dvarw / dmat_par[0][0]
                    else:
                        dvarw = 0.0

                    dmat_par[0][0] = 1.0
                    davew = dvec[0]
                    dvarw = dvarw - (davew ** 2)

                    if dvarw <= 0.0:
                        dvarw = 10**-12

                    for n1 in range(1, ndim + 1):
                        for n2 in range(0, n1):
                            dmat_par[n1][n2] = dmat_par[n2][n1]


                    dmat_par = self.matrix_inv(dmat_par)

                    for n1 in range(0, ndim + 1):
                        for n2 in range(0, ndim + 1):
                            dcoef[n1] = dcoef[n1] + \
                                        dmat_par[n1][n2] * dvec[n2]

                        dpower = dpower + (dcoef[n1] * dvec[n1])

                    dpower = dpower - (davew ** 2)
                    dpowz = (dneff - 3.0) * dpower / (dvarw - dpower) / 2.0
                    dpower = (dneff - 1.0) * dpower / dvarw / 2.0
                    damp = math.sqrt(dcoef[1] * dcoef[1] + \
                                     dcoef[2] * dcoef[2])

                else:

                    dpowz = 0.0
                    dpower = 0.0
                    damp = 0.0

                    if dneff < (10**(-9)):
                        dneff = 0.0

                if damp < (10**(-9)):
                    damp = 0.0

                if dpower < (10**(-9)):
                    dpower = 0.0

                if dpowz < (10**(-9)):
                    dpowz = 0.0

                # Let's write everything out.

                output1.write("%s \t %s \t %s \t %s \t %s \t %s\n" %
                             (str(dtau),str(dfre),str(dpowz),str(damp),
                             str(dcoef[0]),str(dneff)))

                if dpowz > dmz:
                    dmz = dpowz
                    dmfre = dfre
                    dmamp = damp
                    dmcon = dcoef[0]
                    dmneff = dneff

            #
            if max_periods:
                # writes the max_periods output if specified
                output2.write("%f \t %f \t %f \t %f \t %f \t %f\n" %
                             (dtau, dmfre, dmz,  dmamp, dmcon, dmneff))
            
            if gnuplot_compatible:
                # added so that gnuplot reads out of the box
                output1.write("\n")


# If the script runs as a standalone, below is triggered

if __name__ == '__main__':

    # Parsing the arguments

    description = """
    A Weighted Wavelet Z-Transformation Application for Python.
    Translated by M. Emre Aydin - emre.m.aydin@gmail.com
    http://about.me/emre.aydin
    Available at http://github.com/eaydin

    Input arguments can be read from a file. The file descriptor
        prefix is '@'.
    In order to read argument from a file named args.txt,
        the argument @args.txt should be passed.
    An example for args.txt:

        -f=myinputfile.txt
        -o=theoutputfile.output
        -m
        --freq-step=0.001
        -l=0.001
        -hi=0.01
        -c=0.001
        -p=0

    You can pass arguments from file and commandline at the same time.
    If two same arguments passed by this method, the latter will be
    used. So if you want to override some arguments in a an argument
    file, specify the file first.
    An example usage for our earlier @args.txt is as:

        python wwz.py @args.txt -c=0.0125

    The above command will use the settings in args.txt but will
    use c=0.0125 instead of c=0.001

    Comments and blank lines are NOT allowed in argument files.

    Import this script via Python to use it as a module, rather than
    a standalone script. (import wwz)

    """

    parser = argparse.ArgumentParser(prog='wwz.py',
                 formatter_class=argparse.RawDescriptionHelpFormatter,
                 fromfile_prefix_chars="@", description=description)

    parser.add_argument("-f", "--file", type=argparse.FileType("r"),
                        default=sys.stdin, required=True,
                        help="the Input File, Raw Lightcurve")
    parser.add_argument("-o", "--output", type=argparse.FileType('w'),
                        default=sys.stdout, required=True,
                        help="the Output File Name")
    parser.add_argument("-l", "--freq-low", type=float, required=True,
                        help="the Low Frequency Value")
    parser.add_argument("-hi", "--freq-high", type=float, required=True,
                        help="the High Frequency Value")
    parser.add_argument("-d", "--freq-step", type=float, required=True,
                        help="the dF value, incremental step for Frequency")
    parser.add_argument("-c", "--dcon", type=float, required=True,\
                        help="the C constant for the Window Function")
    parser.add_argument("-g", "--gnuplot-compatible", action="store_true",
                        default=False, help="the Output file is GNUPlot \
                        compatible, which means the tau's will be grouped \
                        so that pm3d can easily map. Default value is \
                        'False'.")
    parser.add_argument("-m", "--max-periods", action="store_true",
                        default=False, help="Creates a secondary \
                        output with the maximum Periods for each single \
                        tau. This can be drawn in 2D. The output filename \
                        is derived from the -o option, added 'max_periods'. \
                        Default value is 'False'.")
    parser.add_argument("-t", "--time-divisions", type=float, default=50.0,
                        help="The Time Divisions value. Templeton assumes \
                        this as 50. VStars from AAVSO leaves this optional \
                        contrary to Templeton, yet it's default value is \
                        also 50.")
    parser.add_argument("--time", action="store_true", default=False,
                        help="Calculate the time of operation in seconds \
                        and print to standard output.")
    parser.add_argument("--no-headers", action="store_true", default=False,
                        help="Doesn't print headers to output files if set. \
                        Default is 'False'.")
    parser.add_argument("-p", "--parallel", help="Created threads to speed \
                        up the process. Default value is '1', which means \
                        single thread. '0' means number of detected CPUs,\
                        can be overridden.", type=int, default=1)

    args = parser.parse_args()

    # Check if user asks for time calculation
    if args.time:
        starttime = datetime.now()

    # Get the process number
    cpu=args.parallel

    if cpu != 1:
        try:
            import multiprocessing
        except (ImportError, NotImplementedError):
            print("Multiprocessing not available, using single CPU thread.")
            cpu = 1

    if cpu == 0:
        # detect CPU cores
        try:
            import multiprocessing
            cpu = multiprocessing.cpu_count()
        except (ImportError, NotImplementedError):
            print("Multiprocessing not available, using single CPU thread.")
            cpu = 1

    # Multiprocessing begins        
    if cpu != 1:

        s=WWZPAR(args.file, args.output, args.freq_low, args.freq_high, \
                args.freq_step, args.dcon, args.time_divisions, \
                args.max_periods, args.gnuplot_compatible)

        s.readfile()
        s.maketau()
        s.makefreq()

        # the list of output filenames
        output_par_names = []
        for i in range(1,cpu+1):
            output_par_names.append('wwz.par.proc.%i' % i)

        ntau_par1 = s.ntau/cpu

        # the last ntau_par (turns out this is not necessary)
        ntau_par2 = s.ntau - (ntau_par1 * (cpu-1) )

        thread_list = []

        for i in range(1,cpu+1):
            if i != cpu:
                t = multiprocessing.Process(target=s.wwt, \
                args=(output_par_names[i-1],ntau_par1*(i-1),ntau_par1*i))
            else:
                t = multiprocessing.Process(target=s.wwt, \
                args=(output_par_names[i-1],ntau_par1*(i-1),s.ntau))
            thread_list.append(t)

        for thread in thread_list:
            thread.start()
        for thread in thread_list:
            thread.join()

        # join the outputs and delete remaining files

        if args.max_periods:
            output_max_periods = open(args.output.name + ".max_periods", "w")

        # Write the headers if asked
        if not args.no_headers:
            args.output.write("#%9s %10s %10s %10s %10s %10s\n" % \
                         ("TAU","FREQ","WWZ","AMP","COEF","NEFF"))
            if args.max_periods:
                output_max_periods.write("#%9s %10s %10s %10s %10s %10s\n" % \
                         ("TAU","FREQ","WWZ","AMP","COEF","NEFF"))

        for i in range(1,cpu+1):
            read_par = open(output_par_names[i-1],"r")
            for line in read_par:
                args.output.write(line)
            read_par.close()
            os.remove(output_par_names[i-1])

            if args.max_periods:

                read_par_max = open(output_par_names[i-1] + \
                                    ".max_periods.par","r")
                for line in read_par_max:
                    output_max_periods.write(line)
                read_par_max.close()
                os.remove(output_par_names[i-1] + ".max_periods.par")



    # for single threadding computing
    if cpu == 1:

        # Run the main class and its subroutines

        s=WWZ()

        time_data, magnitude_data = s.readfile(args.file)
        wwz_output = s.wwt(time_data, magnitude_data, args.freq_low, \
                           args.freq_high, args.freq_step, args.dcon, \
                           args.time_divisions)
        if args.gnuplot_compatible:
            send_ntau = len(wwz_output) / \
                  int(((args.freq_high - args.freq_low) / args.freq_step) + 1)

            s.writegnu(wwz_output, args.output, args.no_headers, \
                       args.max_periods, send_ntau)
        else:
            s.writefile(wwz_output, args.output, \
                        args.no_headers, args.max_periods)


    args.output.close()
    if args.max_periods == True and cpu != 1:
        output_max_periods.close()

   # Print calculated time
    if args.time:
        endtime = datetime.now()
        print("Run Time: {0} seconds".format((endtime - starttime).seconds))
