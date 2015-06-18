#-------------------------------------------------------------------------------
# Name:        SAR image simulator
# Purpose:
#
# Author:      vfomin
#
# Created:     08/06/2015
# Copyright:   (c) vfomin 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------


# Numpy, Matplotplib
import numpy as np
import matplotlib.pyplot as plt

# Project
import Tools.GlobalConsts as GC

def logPrint(message):
    print message


LINE_COUNT_LIMIT = 1024*6

class SARSim:

    def __init__(self, platform, verbose=False):
        """

        SAR image simulator
        Geometry :
            - Azimuth direction is along OX axis
            - Range direction is along OY axis

        """
        self._platform = platform
        self._xArray = None
        self.targets = None
        self.verbose = verbose

        print "SAR Platform: ", self._platform.Name

        self._checkSensorParams();


    def setTargets(self, targets):
        """
        Method to set targets
        targets = [[x0,y0,vx0,vy0], [x1,y1,vx1,vy1], ...]

        """
        self.targets = [None]*len(targets)
        for i, target in enumerate(targets):
            assert len(target) == 2 or len(target) == 4, logPrint("A target should be [x0,y0] or [x0,y0,vx0,vy0]")
            self.targets[i] = target if len(target) == 4 else [target[0], target[1], 0.0, 0.0]

    def setParams(self, params, params2=[None, None]):
        """
        Method to set acquisition parameters:
        params = [[xmin, xmax], [theta_min, theta_max], [sceneX, sceneY]]
        params2 = [azimuthSamplingRate (can be None), rangeSamplingRate (can be None)]

        xmin, xmax - azimuth min/max position during the signal acquisition
        theta_min, theta_max - min/max incidence angles (in degrees)
        sceneX, sceneY - zone center position (used for RCMC)
        """
        assert len(params) == 3, logPrint("Aquisition parameters should be of the form [[xmin, xmax], [theta_min, theta_max], [sceneX, sceneY]]")

        self._params = params
        if params2[0] is None:
            dx = self._platform.Vsat / (1.0 * self._platform.PRF)
        else:
            dx = params2[0]
##        dx = self._platform.Vsat / self._platform.PRF

        Nx = int( 1 + (self._params[0][1] - self._params[0][0]) / dx )

        assert Nx > 100, logPrint("Number of azimuth samples is too small. Take a larger azimuth zone ")
        assert Nx < LINE_COUNT_LIMIT, logPrint("Number of azimuth samples, Nx= " + str(Nx) + ", is too large. Take a smaller azimuth zone ")
        self._xArray = np.arange(self._params[0][0], self._params[0][1], dx)
        self._dx = dx

        self._tArray, self._dt = self._computeRangeTimeArray(params2[1])

        print "Number of azimuth samples :", len(self._xArray)
        print "Number of range samples :", len(self._tArray)


    def computeRawData(self):
        """
        Method to compute and get raw data response

        Returns (raw data 2D array, azimuth positions, range time array)
        """
        assert self.targets is not None, logPrint("Targets should be set")
        assert self._xArray is not None and self._tArray is not None, logPrint("Acquisition parameters should be set")
        assert self._dt is not None, logPrint("Error : Simulator badly configured, dt is None")
        assert self._dx is not None, logPrint("Error : Simulator badly configured, dx is None")

        self._rawData = np.zeros((len(self._xArray), len(self._tArray)),dtype=complex)
        deltaX = np.abs(self._params[0][1] - self._params[0][0])
        self._xFootprint = self._computeXFootprint()

        if self.verbose:
            print "XBand: ", deltaX
            print "Azimuth Footprint: ", self._xFootprint

        for target in self.targets:
            dist = lambda x : distance(x / self._platform.Vsat, self._platform.Vsat,
                                self._platform.H, target[0], target[1], target[2], target[3])
            Wa = lambda x : approxWa(x,target[0], self._xFootprint)

            for index, t in enumerate(self._tArray):
                self._rawData[:, index] += received2DSignal(t, self._xArray[:], dist,
                                            self._platform.Tp, self._platform.Freq, self._platform.K, Wa)

        return self._rawData, self._xArray, self._tArray


    def compressRawData2(self, showProgress=True):
        """
        Omega-K compressing algorithm
        1) FFT 2D
        2) Reference Function Multiplication
        3) Stolt interpolation
        4) IFFT 2D
        Returns (compressed 2D array, azimuth positions, range time)
        """
        rawData = self._rawData
        assert rawData is not None, logPrint("Raw data should be computed ( using computeRawData() ) or provided as argument")
        assert self._dt is not None, logPrint("Error : Simulator badly configured, dt is None")
        assert self._dx is not None, logPrint("Error : Simulator badly configured, dx is None")
        assert self._xArray is not None, logPrint("Error : Simulator badly configured, xArray is None")
        assert self._tArray is not None, logPrint("Error : Simulator badly configured, tArray is None")

        Tp = self._platform.Tp
        Vsat = self._platform.Vsat
        K = self._platform.K
        PRF = self._platform.PRF
        H = self._platform.H
        Lambda = self._platform.Lambda
        f0 = self._platform.Freq

        # 1) 2D freq domain
        if showProgress:
            logPrint(" - 2D FFT ...")

        fRangeArray = GC.freq(rawData.shape[1], 1.0/self._dt)
        fAzimuthArray = GC.freq(rawData.shape[0], Vsat/self._dx)

        procData = rawData
        procData = GC.fft2(procData)

        if self.verbose:
            showResponse(procData, None, '2D FFT')


        # 2) RFM
        if showProgress:
            logPrint(" -- RMF ...")

        # Reference distance
        tZero = self._tArray[len(self._tArray)/2]
        etaZero = self._xArray[len(self._xArray)/2] / Vsat
        Rref = 0.5 * GC.C * tZero
        Vref = Vsat
        phaseRef = generateRefPhase(fAzimuthArray, fRangeArray, f0, Rref, Vref, K, Tp, tZero, etaZero)

        procData = procData * phaseRef
        if self.verbose:
##            showReImAbsPhaseResponse(phaseRef, None, 'Reference Phase')
            showResponse(procData, None, 'RFM')

        # 3) Stolt mapping
        if showProgress:
            logPrint(" --- Stolt mapping ...")

        # For azimuth line :
        for i, f in enumerate(fAzimuthArray):
            #
            # Here we need analytical inverse of Stolt mapping :
            # f' + f0 = sqrt( (f0 + f)^2 - (c * fEta / (2*Vsat))^2 )
            #
            # f + f0 = + sqrt( (f0 + f')^2 + (c * fEta / (2*Vsat))^2 )
            #
            nfRangeArray = -f0 + np.sqrt( (f0 + fRangeArray)**2 + ( 0.5 * GC.C * f / Vsat)**2 )

            yRe = np.real(procData[i,:])
            yIm = np.imag(procData[i,:])
            yRe = np.interp(nfRangeArray, fRangeArray, yRe)
            yIm = np.interp(nfRangeArray, fRangeArray, yIm)
            procData[i,:] = yRe + 1j * yIm

        if self.verbose:
            showReImAbsPhaseResponse(procData, None, 'Stolt mapping')

        # 4) 2D IFFT
        if showProgress:
            logPrint(" ---- 2D IFFT ...")

        procData = GC.ifft2(procData)

        if self.verbose:
            showResponse(procData, None, '2D IFFT')

        return procData, self._xArray, self._tArray


    def compressRawData(self, showProgress=True):
        """
        Apply RDA on raw data
        1) Range compression
        2) Range Cell Migration Correction
        3) Azimuth compression

        Returns (compressed 2D array, azimuth positions, range time)
        """
        rawData = self._rawData
        assert rawData is not None, logPrint("Raw data should be computed ( using computeRawData() ) or provided as argument")
        assert self._dt is not None, logPrint("Error : Simulator badly configured, dt is None")
        assert self._dx is not None, logPrint("Error : Simulator badly configured, dx is None")
        assert self._xArray is not None, logPrint("Error : Simulator badly configured, xArray is None")
        assert self._tArray is not None, logPrint("Error : Simulator badly configured, tArray is None")

        Tp = self._platform.Tp
        Vsat = self._platform.Vsat
        K = self._platform.K
        PRF = self._platform.PRF
        H = self._platform.H
        Lambda = self._platform.Lambda

        # 1) Range compression :
        if showProgress:
            logPrint(" - Range compression ...")

        fRangeArray = GC.freq(rawData.shape[1], 1.0/self._dt)
        HMF = GC.rect(fRangeArray / (np.abs(K)*Tp)) * np.exp(1j*np.pi*fRangeArray**2 / K)
        # Multiply by a phase due to different range zero time
##        tZero = self._tArray[len(self._tArray)/2]
##        HMF = HMF * np.exp( -1j* 2.0 * np.pi * fRangeArray * tZero)
        procData = rangeCompression(rawData, HMF)
        HMF = None

        if self.verbose:
            showResponse(procData, None, 'Range compression')

        # 2) Azimuth FFT :
        if showProgress:
            logPrint(" -- Azimuth FFT ...")
        procDataAF = np.zeros(procData.shape, dtype=complex)
        for i in range(procDataAF.shape[1]):
            # get range fixed column
            procDataAF[:,i] = GC.fft(procData[:,i])
        procData = procDataAF
        procDataAF = None


        fArray = GC.freq(procData.shape[0], self._platform.Vsat/self._dx)

        if self.verbose:
##            showResponse(procData, [self._tArray[0], self._tArray[-1], fArray[0], fArray[-1]], 'Azimuth FFT')
            showResponse(procData, None, 'Azimuth FFT')

        # 3) Range migration :
        if showProgress:
            logPrint(" --- Range migration ...")
        Dfreq = np.sqrt(1.0 - (0.5 * Lambda * fArray / Vsat)**2)
        procData = rangeMigration2(procData, self._tArray, Dfreq)

        if self.verbose:
##            showResponse(procData, [self._tArray[0], self._tArray[-1], fArray[0], fArray[-1]], 'RCMC')
            showResponse(procData, None, 'RCMC')


        # 4) Azimuth compression :
        if showProgress:
            logPrint(" ---- Azimuth compression ...")

        etaZero = self._xArray[len(self._xArray)/2] / Vsat
        for i in range(procData.shape[1]):
            # Use H(f) = exp( 1j * 4 * pi * R0 * D(f) * f0 / c )
            #
            HMF = np.exp( 2*np.pi*1j * GC.C * self._tArray[i] / Lambda * Dfreq)
            # Multiply by a phase due to different azimuth zero time
            HMF = HMF * np.exp( -1j* 2.0 * np.pi * fArray * etaZero)
            procData[:,i] = procData[:,i] * HMF


        # 5) Azimuth IFFT
        if showProgress:
            logPrint(" ----- Azimuth IFFT ...")
        procDataAF = np.zeros(procData.shape, dtype=complex)
        for i in range(procDataAF.shape[1]):
            # get range fixed column
            procDataAF[:,i] = GC.ifft(procData[:,i])
        procData = procDataAF
        procDataAF = None

        if self.verbose:
            showResponse(procData, [self._tArray[0], self._tArray[-1], fArray[0], fArray[-1]], 'Compressed image', False)


        return procData, self._xArray, self._tArray




    def _computeXFootprint(self):
        thetaX = 0.886 * self._platform.Lambda / self._platform.La
        angleC = (self._params[1][1] + self._params[1][0])*0.5
        ycenter = self._platform.H * np.tan(angleC * np.pi / 180.0)
        rc = np.sqrt(self._platform.H**2 + ycenter**2)
        return 2.0 * np.tan(thetaX *0.5) * rc

    def _computeRangeTimeArray(self, rangeSamplingRate=None):

        if rangeSamplingRate is None:
            B = self._platform.Tp*np.abs(self._platform.K)
            dt = 1.0/(1.4*B)
        else:
            dt = 1.0/rangeSamplingRate

        ymin = self._platform.H * np.tan(self._params[1][0] * np.pi / 180.0)
        ymax = self._platform.H * np.tan(self._params[1][1] * np.pi / 180.0)
        rmin = np.sqrt(self._platform.H**2 + ymin**2)
        deltaX = np.abs(self._params[0][1] - self._params[0][0])
        rmax = np.sqrt(self._platform.H**2 + deltaX**2 + ymax**2)
        if self.verbose:
            print "rmin, rmax", rmin, rmax

        tmin = -0.50*self._platform.Tp + 2.0 * rmin / GC.C
        tmax = 0.50*self._platform.Tp + 2.0 * rmax / GC.C
##        tmin = 2.0 * rmin / GC.C
##        tmax = 2.0 * rmax / GC.C

        if self.verbose:
            print "tmin, tmax", tmin, tmax

        Nt = int( 1 + (tmax - tmin) / dt )
        assert Nt > 100, logPrint("Number of range samples, Nt=" + str(Nt) + ", is too small.")
        assert Nt < LINE_COUNT_LIMIT, logPrint("Number of range samples, Nt=" + str(Nt) + ", is too large. Take a smaller acquisition angles")

        tArray = np.arange(tmin, tmax, dt)
        if (len(tArray) % 2 == 1):
            tArray = np.resize(tArray, (len(tArray)-1))
        return tArray, dt

    def _checkSensorParams(self):
        PRF = self._platform.PRF
        Lambda = self._platform.Lambda
        Vsat = self._platform.Vsat
        assert self._platform.PRF * 0.5 < 2.0 / Lambda * Vsat, logPrint("Frequency azimuth max (PRF/2) is larger than 2/Lambda * Vsat ")



def generateRefPhase(fAzimuthArray, fRangeArray, f0, Rref, Vref, K, Tp, tZero=0, etaZero=0):
    """
    Phi_ref(f_eta, f_tau) = 4*pi*Rref/c * sqrt( (f0+f_tau)^2 - c^2 * f_eta^2 / (4*Vref^2) )
    - tZero corresponds to the 'zero' range time of the signal range time sampling. tZero = (tmax+tmin)/2
    - etaZero corresponds to the 'zero' azimuth time of the signal azimuth time sampling. etaZero = (etaMax+etaMin)/2

    """
    a = 4.0*np.pi*Rref/GC.C
    fRangeMax = np.abs(K)*Tp
##        fRangeMax = fRangeArray[-1]
    phase = np.zeros((len(fAzimuthArray),len(fRangeArray)), dtype=np.complex)
    for i, f in enumerate(fAzimuthArray):
        b = np.sqrt( (f0 + fRangeArray)**2 - ( 0.5 * GC.C * f / Vref)**2 )
        c = np.pi * fRangeArray**2 / K
        d =  2.0*np.pi* tZero * fRangeArray
        d2 = 2.0*np.pi* etaZero * f
        phase[i,:] = GC.rect(fRangeArray / fRangeMax) * np.exp(1j * a * b[:] + 1j*c[:] - 1j*d[:] - 1j*d2)
    return phase



def rangeCompression(Srx, HMF):
    """
    Range compression with freq-domain filter H(f)
    """
    assert Srx.shape[1] == len(HMF), logPrint("Range dimension of the Srx is not equal to length of HMF")

    SrxMF=np.zeros(Srx.shape,dtype=complex)
    for i in range(Srx.shape[0]):
        data = Srx[i,:]
        DataMF = GC.fft(data) * HMF
        SrxMF[i,:] = GC.ifft(DataMF)

    return SrxMF



def compression(Srx, h, dim0):
    """
    Matched filter alogorith with filter 'h' is applied on
    the image 'Srx' on the dimension 'dim'

    dim is 'range' or 'azimuth'

    """

    dim = 1 if dim0 is 'range' else 0 if dim0 is 'azimuth' else -1
    assert dim >= 0, logPrint("Parameter dim should be 'range' or 'azimuth'")
    # inverse of dim
    idim = np.abs(dim-1)
    # pad received and emitted signals on common time stemp
    # to obtain common fft
    n1 = len(h)
    n2 = Srx.shape[dim]

    Nsa = n1 + n2 - 1
##    Nsa2 = int(np.power(2.0, np.ceil(np.log2( Nsa ))))
    Nsa2 = Nsa if Nsa % 2 == 0 else Nsa+1


##    shape = (Srx.shape[0], Nsa2) if dim == 1 else (Nsa2, Srx.shape[1])
##    SrxMF=np.zeros(shape,dtype=complex)
    SrxMF=np.zeros(Srx.shape,dtype=complex)
    hMF = np.concatenate((np.zeros(np.floor((Nsa2-n1)*0.5)), h, np.zeros(np.ceil((Nsa2-n1)*0.5))))
    HMF = GC.fft(hMF)

    m1 = np.floor((Nsa2-n2)*0.5)
    m2 = np.ceil((Nsa2-n2)*0.5)
    for i in range(Srx.shape[idim]):
        data = Srx[i,:] if dim == 1 else Srx[:,i]
        dataMF = np.concatenate((np.zeros(m1), data, np.zeros(m2)))
        DataMF = GC.fft(dataMF)
        DataMF = DataMF * HMF
        if dim == 1:
            SrxMF[i,:] = GC.ifft(DataMF)[n1/2:n1/2+n2]
        else:
            SrxMF[:,i] = GC.ifft(DataMF)[n1/2:n1/2+n2]

    return SrxMF


def deltaR(x,xc,dist):
    """
    delta R for a stationary ground target zone
    deltaR = R(x) - Rc,
    where R(x) is the distance when sensor is at azimuth point x
    Rc=(xc,yc) is the smallest distance to the center of the zone

    - x is the numpy array of azimuth positions
    - xc is the azimuth position of the targeted zone
    - dist is the distance function

    """
    return dist(x) - dist(xc)

def rangeMigration(Srx, dRI):
    """
    RCMC algorithm
    """
    Srm = np.zeros(Srx.shape, dtype=complex)
    rcm_max = dRI.max()

    # Loop on rows
    ind1=np.arange(int(Srm.shape[1] - rcm_max),dtype=int)
    for i in range(Srm.shape[0]):
        # Loop on row elements:
        ind2 = ind1+dRI[i]
        Srm[i,ind1[:]] = Srx[i,ind2]

    return Srm

def rangeMigration2(Srx, tArray, Dfreq):
    """
    RCMC algorithm with interpolation
    """
    assert Srx.shape[0] == len(Dfreq), logPrint("Error : Srx.shape[0] != len(Dfreq)")
    assert Srx.shape[1] == len(tArray), logPrint("Error : Srx.shape[1] != len(tArray)")

    Srm = np.zeros(Srx.shape, dtype=complex)

    for i in range(Srm.shape[0]):
        tArray2 = tArray / Dfreq[i]
        # Loop on row elements:
        yRe = np.real(Srx[i,:])
        yIm = np.imag(Srx[i,:])

        yRe = np.interp(tArray2, tArray, yRe)
        yIm = np.interp(tArray2, tArray, yIm)

        Srm[i,:] = yRe + 1j * yIm

    return Srm


def distance(etaArray, Vsat, H, x0, y0, vx0=0.0, vy0=0.0, ax0=0.0, ay0=0.0):
    """
    Range distance between sensor and a single ground target
    R(etaArray)^2 = H^2 + (y0 + vy0 * etaArray + 0.5*ay0 * etaArray^2)^2
                        + (x0 + vx0 * etaArray + 0.5*ax0 * etaArray^2 - Vsat * etaArray)^2
    x0, y0 is target's position
    vx0, vy0 is target's velocity
    ax0, ay0 is target's acceleration
    Vsat is sensor's velocity
    etaArray is the azimuth time
    """
    y2 = y0 + vy0 * etaArray + 0.5*ay0*etaArray**2
    y2 = y2**2
    x2 = x0 + vx0 * etaArray + 0.5*ax0*etaArray**2 - Vsat * etaArray
    x2 = x2**2
    return np.sqrt(H**2 + y2 + x2)

def received2DSignal(t, x, dist, tau_p, f_c, K_r, Wa, sigma=1) :
    """
    SAR received signal at range time t and azimuth position x
    S_{out}(t,x) = rect( (t-2*d(x)/c) / tau_p) * Wa(x) * sigma
                    * exp( -4*pi*1j* f_c * d(x)/c + 1j*pi*K_r*(t-2*d(x)/c)^2 )

    """
    tt = t - 2.0*dist(x) / GC.C
    phase = -4.0*f_c*dist(x) / GC.C + K_r*(tt**2)
    return GC.rect( tt/tau_p ) * Wa(x) * sigma * np.exp(1j*np.pi*phase)


def approxWa(x,x_c,deltaX):
    """
    Antenna Beam pattern
    Approximation of a uniform aperture (see SARAS: A synthetic Aperture Radar (SAR) Raw Signal Simulator, G.Franceschetti)
    Wa(x) = sinc^2 ( (x - x_c)/deltaX ) where
        deltaX - azimuth swath corresponding to -3.9dB points
    """
    wa = np.sinc( (x-x_c) / deltaX )
    return wa**2

def showReImAbsPhaseResponse(Srx, extent=None, title="Figure", keepAspectRatio=True, ):
    f = plt.figure()
    f.suptitle(title)
    plt.subplot(221)
    plt.title("SAR response, module")
    plt.xlabel("Range")
    plt.ylabel("Azimuth")
    plt.imshow(np.absolute(Srx), extent=extent, aspect='auto' if not keepAspectRatio else None, interpolation='none')
    plt.colorbar()
    plt.subplot(222)
    plt.title("SAR response, phase")
    plt.xlabel("Range")
    plt.ylabel("Azimuth")
    plt.imshow(np.angle(Srx), extent=extent, aspect='auto' if not keepAspectRatio else None, interpolation='none')
    plt.colorbar()
    plt.subplot(223)
    plt.title("SAR response, Re")
    plt.xlabel("Range")
    plt.ylabel("Azimuth")
    plt.imshow(np.real(Srx), extent=extent, aspect='auto' if not keepAspectRatio else None, interpolation='none')
    plt.colorbar()
    plt.subplot(224)
    plt.title("SAR response, Im")
    plt.xlabel("Range")
    plt.ylabel("Azimuth")
    plt.imshow(np.imag(Srx), extent=extent, aspect='auto' if not keepAspectRatio else None, interpolation='none')
    plt.colorbar()
    plt.show()

def showGrayResponse(Srx, extent=None, title="Figure", keepAspectRatio=True):
    f = plt.figure()
    f.suptitle(title)
    plt.title("SAR response, module")
    plt.xlabel("Range")
    plt.ylabel("Azimuth")
    plt.imshow(np.absolute(Srx), extent=extent, aspect='auto' if not keepAspectRatio else None, interpolation='none')
    plt.colorbar()
    plt.set_cmap('gray')
    plt.show()
    plt.set_cmap('jet')

def showResponse(Srx, extent=None, title="Figure", keepAspectRatio=True):
    f = plt.figure()
    f.suptitle(title)
    plt.title("SAR response, module")
    plt.xlabel("Range")
    plt.ylabel("Azimuth")
    plt.imshow(np.absolute(Srx), extent=extent, aspect='auto' if not keepAspectRatio else None, interpolation='none')
    plt.colorbar()
    plt.show()


def warpFromRangeTimeToRegularOYGrid(Srx, tArray, H):
    """
    Method to warp data from range time grid into a regular y grid
    y(t) = sqrt( (1/2 * C * t)^2 - H^2 )

    """

    yArray = np.sqrt( (0.5*GC.C*tArray)**2 - H**2 )

##    dy = (yArray[-1] - yArray[0])/(len(yArray)-1)
##    yArray2 = np.arange(yArray[0],yArray[-1],dy)
    yArray2 = np.linspace(yArray[0],yArray[-1],len(yArray))


    Sout = np.zeros(Srx.shape, dtype=np.complex)
    for i in range(Srx.shape[0]):
        yRe = np.real(Srx[i,:])
        yIm = np.imag(Srx[i,:])
        yRe = np.interp(yArray2, yArray, yRe)
        yIm = np.interp(yArray2, yArray, yIm)
        Sout[i,:] = yRe + 1j*yIm

    return Sout, yArray2







if __name__ == '__main__':

    import Sensors.GlobalVariablesCSK_HI as Sensor
    # Test
    xc = 0.0
    yc = Sensor.H * np.tan(55.0 * np.pi / 180.0)

    # Absolute reference distance
    Rref = np.sqrt(Sensor.H**2 + yc**2)
    # Relative reference distance to Rmin
    ymin = Sensor.H * np.tan(54.995 * np.pi / 180.0)
    Rref2 = Rref - np.sqrt(Sensor.H**2 + ymin**2)
    Vref = Sensor.Vsat

    dt = 1.0/(1.4*np.abs(Sensor.K) * Sensor.Tp)
    dx = Vref/Sensor.PRF

    fRangeArray = GC.freq(2*1024, 1.0/dt)
    fAzimuthArray = GC.freq(2*1024, Vref/dx)

    phaseRef = generateRefPhase(fAzimuthArray, fRangeArray, Sensor.Freq, Rref2, Vref, Sensor.K, Sensor.Tp)

    Sref = GC.ifft2(phaseRef)
    showResponse(Sref)


    exit()

##    filename = "C:\\VFomin_folder\\PISE_project\\SAR-GMTI\\Data\\CSK_RAW_0340_B001_2689_15971_3028_3027.npy"
##    pt = [2689,15971]

##    filename = "C:\\VFomin_folder\\PISE_project\\SAR-GMTI\\Data\\CSK_RAW_0340_B001_10255_8742_3195_3195.npy"
##    pt = [10255,8742]

##    filename = "C:\\VFomin_folder\\PISE_project\\SAR-GMTI\\Data\\CSK_RAW_0340_B001_19165_15130_3532_3364.npy"
##    pt = [19165,15130]

##    data = np.load(filename)
##    print data.shape, data.dtype
##    # (x,y,w,h) and x=range, y=azimuth
##    win = [pt[0],pt[1],data.shape[1],data.shape[0]]
##    showResponse(data)
##
##
##    # TEST USING CSK_RAW data
##    import GlobalVariablesCSK_HI as Sensor
##    Sensor.K = -2.1661376953124998e12
##    Sensor.H = 624840.9574418962
##    SamplingRate = 1.125e8
##    time_bias = 14612
##    ECEFSatelliteVelocity = "-3340.3819748,-106.1221533,-6861.7438606,-3275.9787599,-63.7058631,-6893.3510286,-3211.1342775,-21.376226,-6924.1548523,-3145.8562251,20.86118,-6954.1516992,-3080.1523606,63.0007979,-6983.3380332,-3014.0305018,105.0370925,-7011.7104149,-2947.4985255,146.9645515,-7039.2655026,-2880.5643662,188.7776861,-7066.0000528,-2813.2360152,230.4710319,-7091.9109206,-2745.5215193,272.0391502,-7116.9950602,-2677.4289799,313.4766287,-7141.2495257,-2608.9665514,354.7780829,-7164.6714711,-2540.1424401,395.9381564,-7187.2581505,-2470.9649028,436.9515219,-7209.0069185,-2401.4422456,477.8128826,-7229.9152297"
##    Sensor.Vsat = Sensor._computeMeanVelocity(ECEFSatelliteVelocity,',')
##    AbsoluteAzimuthFirstTime = 14612.406045371094 - time_bias
##    AbsoluteAzimuthLastTime = 14619.758871022656 - time_bias
##    RangeFirstTimes = 0.005281833666666667
##    azimuthSamplingRate = 1.0 / Sensor.PRF
##    rangeSamplingRate = 1.0/ SamplingRate
##
##    # y=azimuth, x=range
##    winAzimuthFirstTime = AbsoluteAzimuthFirstTime + win[1] * azimuthSamplingRate
##    winAzimuthLastTime = winAzimuthFirstTime + win[3] * azimuthSamplingRate
##    winRangeFirstTime = RangeFirstTimes + win[0] * rangeSamplingRate
##    winRangeLastTime = winRangeFirstTime + win[2] * rangeSamplingRate
##
##    s = SARSim(Sensor,verbose=False)
##    s._rawData = data
##
##    s._dt = 1.0/SamplingRate
##    Nt = int( 1 + (winRangeLastTime - winRangeFirstTime) / s._dt )
##    assert Nt > 100, logPrint("Number of range samples is too small.")
##    assert Nt < 4096, logPrint("Number of range samples is too large. Take a smaller acquisition angles")
##    s._tArray = np.arange(winRangeFirstTime, winRangeLastTime, s._dt)
##    if len(s._tArray) > data.shape[1]:
##        s._tArray = np.resize(s._tArray, (data.shape[1]))
##
##    s._dx = Sensor.Vsat / (1.0 * Sensor.PRF)
##    Nx = int( 1 + Sensor.Vsat * (winAzimuthLastTime - winAzimuthFirstTime) / s._dx )
##    assert Nx > 100, logPrint("Number of azimuth samples is too small. Take a larger azimuth zone ")
##    assert Nx < 4096, logPrint("Number of azimuth samples is too large. Take a smaller azimuth zone ")
##    s._xArray = Sensor.Vsat * np.arange(winAzimuthFirstTime, winAzimuthLastTime, azimuthSamplingRate)
##    if len(s._xArray) > data.shape[0]:
##        s._xArray = np.resize(s._xArray, (data.shape[0]))
##
####    s._xSwath = s._xArray[-1] - s._xArray[0]
##
##    cdata, xArray, tArray = s.compressRawData()
##    showResponse(cdata)
##
##    plt.imshow(np.abs(cdata))
##    plt.gray()
##    plt.show()
##    exit()




    TOL = 1e-10
    # Test Aerial
    import Sensors.GlobalVariablesAerial as Sensor
    # Test
    xc = 0.0
    yc = Sensor.H * np.tan(55.0 * np.pi / 180.0)
    s = SARSim(Sensor,verbose=False)
    # define SAR acquisition parameters:
    params = [[-500, 500], [52.0, 57.0], [xc, yc]]
    s.setParams(params)
    # set targets

    # a) stationary targets :
    setOfTargets = [
        [[xc, yc]],
        [[xc + 123, yc - 213]],
        [[xc, yc], [xc + 123, yc - 213]],
    ]

    S = [None]*len(setOfTargets)
    for i, targets in enumerate(setOfTargets):
        s.setTargets(targets)
        Srx, xArray, tArray = s.computeRawData()
        S[i], xArray, tArray = s.compressRawData()

    # compare results:
    diff = S[0] + S[1] - S[2]
    diff = np.abs(diff)**2
    s = np.sum(diff)
    print s, s.min(), s.max()

    assert s.max() < TOL, logPrint("Processing error: Sum(el) != Sum")


    # b) moving targets :
##    setOfTargets = [
##        [[xc, yc]],
##        [[xc + 123, yc - 213, 0.0, 2.0]],
##        [[xc - 123, yc + 213, 3.0, 0.0]],
##        [[xc + 123, yc + 213, 2.1, 1.7]],
##        [[xc, yc], [xc + 123, yc - 213, 0.0, 2.0], [xc - 123, yc + 213, 3.0, 0.0], [xc + 123, yc + 213, 2.1, 1.7]]
##    ]
##
##    S = [None]*len(setOfTargets)
##    for i, targets in enumerate(setOfTargets):
##        s.setTargets(targets)
##        Srx, xArray, tArray = s.computeRawData()
##        S[i], xArray, tArray = s.compressRawData()
##
##    # compare results:
##    diff = S[0] + S[1] + S[2] + S[3] - S[4]
##    diff = np.abs(diff)**2
##    s = np.sum(diff)
##    print s, s.min(), s.max()
##
##    assert s.max() < TOL, logPrint("Processing error: Sum(el) != Sum")

    exit()


    # Test CSK
    import Sensors.GlobalVariablesCSK_HI as Sensor
    # Test
    xc = 0.0
    yc = Sensor.H * np.tan(51.00005 * np.pi / 180.0)
    s = SARSim(Sensor,verbose=True)
    # define SAR acquisition parameters:
    params = [[-3000, 3000], [51.0, 51.00010], [xc, yc]]
    s.setParams(params)
    # set targets
    targets = [[xc + 100.0, yc + 100.0], [xc - 150.0, yc - 50.0], [xc, yc], [xc + 200.0, yc]]
##    targets = [[xc, yc], [xc + 100.0, yc], [xc - 100.0, yc]]
##    targets = [[xc, yc]]
    s.setTargets(targets)

    Srx, xArray, tArray = s.computeRawData()
    showResponse(Srx, [tArray[0],tArray[-1],xArray[0],xArray[-1]], 'Raw data')

    exit()

    # data compression
    S_final, xArray, tArray = s.compressRawData()
    showResponse(S_final, [tArray[0],tArray[-1],xArray[0],xArray[-1]], 'Compressed image')

    exit()

