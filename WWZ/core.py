import math
import sys
import os
from datetime import datetime
import numpy


class WWZ:

    def __init__(self):
        self.wwz = []

    @staticmethod
    def makefreq(flo, fhi, df):
        # Lines 149 - 181 Fortran
        # Lines 356 - 370 Java

        freq = [flo]
        nfreq = int((fhi - flo)/ df) + 1

        # These lines seem skeptical!
        for i in range(1, nfreq+1):
            freq.append(flo + ((i - 1)* df))
        return freq

    @staticmethod
    def roundtau(darg):
        dex = math.log(darg, 10)
        nex = int(dex)

        darg = darg / math.pow(10, nex)

        if darg >=5:
            darg = 5.0
        elif darg >= 2:
            darg = 2.0
        else:
            darg = 1.0

        darg = darg * math.pow(10, nex)
        return darg

    def maketau(self, time_series, timedivisions):
        # The Maketau section
        # Lines 90 - 122 Fortran

        dtauhi = time_series[-1]
        dtaulo = time_series[0]

        dtspan = dtauhi - dtaulo
        dtstep = self.roundtau(dtspan / timedivisions)
        dtaulo = dtstep * int(dtaulo / dtstep)
        dtauhi = dtstep * int((dtauhi / dtstep) + 0.5)

        tau = []
        dtau = dtaulo
        while dtau <= dtauhi:
            tau.append(dtau)
            dtau = dtau + dtstep
        return tau


    @staticmethod
    def matrix_inversion(input_matrix):
        # Lines 202 - 252 Fortran

        ndim = 2
        dsol = numpy.zeros(shape=(3, 3))

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

    def wwt(self, time_series, magnitude, flo, fhi, df, dcon, timedivisions):

        freq = self.makefreq(flo, fhi, df)
        nfreq = len(freq)
        dmat = numpy.zeros(shape=(3, 3))

        tau = self.maketau(time_series, timedivisions)
        ntau = len(tau)

        dvec = [0, 0, 0]
        dcoef = [0, 0, 0]

        dvarw = 0.0

        twopi = 2.0 * math.pi
        ndim = 2
        itau2 = ntau
        ifreq1 = 1
        ifreq2 = nfreq

        # Creating output arrays
        output = numpy.empty((ntau*(nfreq-1), 6))
        numdat = len(time_series)
        index = 0

        for itau in range(0, itau2):
            nstart = 1
            dtau = tau[itau]

            dmz = -1.0  # less than the smallest WWZ
            for ifreq in range(ifreq1, ifreq2):
                dfre = freq[ifreq]
                domega = dfre * twopi

                # THESE LINES CAN BE BETTER?
                for i in range(0, ndim + 1):
                    dvec[i] = 0.0
                    for j in range(0, ndim + 1):
                        dmat[i][j] = 0.0
                dweight2 = 0.0
                # --------------------------


                for idat in range(nstart, numdat):
                    dz = domega * (time_series[idat] - dtau)
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

                # THIS CAN BE FASTER
                for n1 in range(0, ndim + 1):
                    dcoef[n1] = 0.0
                # ----------------------

                if dweight2 > 0:
                    dneff = (dmat[0][0] * dmat[0][0]) / dweight2
                else:
                    dneff = 0.0

                if dneff > 3.0:

                    # THIS CAN BE FASTER
                    for n1 in range(0, ndim + 1):
                        dvec[n1] = dvec[n1] / dmat[0][0]

                        for n2 in range(1, ndim + 1):
                            dmat[n1][n2] = dmat[n1][n2] / dmat[0][0]
                    # -------------
                    if dmat[0][0] > 0.005:
                        dvarw = dvarw / dmat[0][0]
                    else:
                        dvarw = 0.0

                    dmat[0][0] = 1.0
                    davew = dvec[0]
                    dvarw = dvarw - (davew ** 2)

                    if (dvarw < 0):
                        dvarw = 10**-12

                    for n1 in range(1, ndim + 1):
                        for n2 in range(0, n1):
                            dmat[n1][n2] = dmat[n2][n1]

                    dmat = self.matrix_inversion(dmat)

                    for n1 in range(0, ndim + 1):
                        for n2 in range(0, ndim + 1):
                            dcoef[n1] = dcoef[n1] + dmat[n1][n2] * dvec[n2]
                        dpower = dpower + (dcoef[n1] * dvec[n1])

                    dpower = dpower - (davew ** 2)
                    dpowz = (dneff - 3.0) * dpower / (dvarw - dpower) / 2.0
                    dpower = (dneff - 1.0) * dpower / dvarw / 2.0
                    damp = math.sqrt(dcoef[1] * dcoef[1] + dcoef[2] * dcoef[2])
                else:
                    dpowz = 0.0
                    dpower = 0.0
                    damp = 0.0

                if dneff < (10 ** (-9)):
                    dneff = 0.0

                if damp < (10 ** (-9)):
                    damp = 0.0

                if dpower < (10 ** (-9)):
                    dpower = 0.0

                if dpowz < (10 ** (-9)):
                    dpowz = 0.0

                output[index] = [dtau, dfre, dpowz, damp, dcoef[0], dneff]
                index += 1
                if dpowz > dmz:
                    dmz = dpowz
                    dmfre = dfre
                    dmamp = damp
                    dmcon = dcoef[0]
                    dmneff = dneff

        self.wwz = output

    def write_file(self, fp, headers=False):
        numpy.set_printoptions(precision=5)
        numpy.set_printoptions(suppress=True)
        numpy.set_printoptions(threshold='nan')

        if headers:
            numpy.savetxt(fp, self.wwz, delimiter='\t',
                          fmt="%10.4f", comments='#',
                          header="{:9s} {:10s} {:10s} {:10s} {:10s} {:10s}"
                          .format("TAU", "FREQ", "WWZ", "AMP", "COEF", "NEFF"))
        else:
            numpy.savetxt(fp, self.wwz, delimiter='\t', fmt="%10.4f")
