#!/usr/bin/env python
"""
This plot displays the audio waveform, spectrum, and spectrogram from the 
microphone.

Based on updating_plot.py
"""
import sys

# Major library imports
try:
    import pyaudio
except ImportError:
    sys.exit('You need pyaudio installed to run this demo.')

from numpy import zeros, linspace, short, fromstring, hstack, transpose
from scipy import fft
from scipy.signal import blackmanharris, butter, filtfilt
import numpy as np

# Enthought library imports
from chaco.default_colormaps import hot
from enable.api import Component, ComponentEditor
from traits.api import HasTraits, Instance
from traitsui.api import Item, Group, View, Handler
from pyface.timer.api import Timer

# Chaco imports
from chaco.api import Plot, ArrayPlotData, HPlotContainer


SAMPLING_RATE = 11025
FRAMES_PER_PLOT = 100
OCTAVES = 3 
TOP_MIDINUMBER = 85
BINS = 12

def sqrt_blackmanharris(N):
    return np.sqrt(blackmanharris(N))

def nextpow2(n):
    p2 = 2
    while p2 < n:
        p2 *= 2
    return p2

class CQTKernel(object):
   
    def __init__(self, max_midinumber, bins, fs,
                 q = 6.4,
                 atomHopFactor = 0.25,
                 thresh = 0.0005, 
                 winFunc = sqrt_blackmanharris,
                 ):
        """
        Generate the CQT kernel for one octave
        """
        fmax = 2**((max_midinumber - 69) / 12.) * 440
        if fmax >= fs/2.:
            raise ValueError("fmax (%s) is too big for fs (%s)"
                             %(str(fmax), str(fs)))
        F = fmax * 2 ** ((np.arange(bins) - bins + 1) / 12.)  
        fmin = F[0]
        Q = q / (2 ** (1./bins) - 1)
        N = np.round(Q * fs / F)                             
        
        atomHOP = nextpow2(N[-1] * atomHopFactor)/2
        first_center = atomHOP * np.ceil((N[0] / 2.) / atomHOP)
        FFTLen = nextpow2(first_center + np.ceil(N[0] / 2.))
        atomNr = np.floor((FFTLen - np.ceil(N[0] / 2.) -
                              first_center) / atomHOP) + 1 
        last_center = first_center + (atomNr - 1.) * atomHOP
        fftHOP = (last_center + atomHOP) - first_center
        
        sparKernel = np.zeros([bins * atomNr, FFTLen],
                              dtype=np.complex)
        
        
        # Compute kernel
        for k in np.arange(bins): 
            
            Nk = N[k]
            window = winFunc(Nk)
            
            fk = F[k]
            tempKernelBin = (window *1. / Nk) * np.exp(2 * np.pi * 1j * fk * np.arange(Nk) / fs)
            atomOffset = first_center - np.ceil(Nk / 2.)
            tempKernel = np.zeros(FFTLen, dtype=np.complex)
            for i in np.arange(atomNr): 
                shift = atomOffset + i * atomHOP
                tempKernel[shift:(Nk+shift)] = tempKernelBin 
                sparKernel[int(i + k * atomNr)] = np.fft.fft(tempKernel)
                tempKernel = np.zeros(FFTLen, dtype=np.complex)  
        
        sparKernel[np.abs(sparKernel)<=thresh] = 0
        
        wx1 = np.argmax(sparKernel[0,:])
        wx2 = np.argmax(sparKernel[-1,:])
        wK = sparKernel[:,wx1:wx2]
        wK = np.diag(np.dot(wK.T, np.conjugate(wK)))
        wK = wK[int(np.round(1./q)):\
                int(len(wK) - np.round(1./q) - 1)]
        weight = 1. / np.mean(np.abs(wK))
        weight *= (fftHOP * 1. / FFTLen)
        weight = np.sqrt(weight)
        sparKernel *= weight

        K = np.conjugate(sparKernel)

        self.B, self.A = butter(N=6, Wn=0.5, btype='low')
        self.FFTLen = FFTLen
        self.bins = bins
        self.atomNr = atomNr
        self.K = K
       
    def __str__(self):
        description = "CQT Kernel structure, containing:\n"
        for k, v in self.__dict__.items():
            description += str(k) + ': ' + str(v) + '\n'
        return description

    def cqt_oneOct(self,x):
        assert x.size >= self.FFTLen
        X = np.fft.fft(x[-self.FFTLen:], n=self.FFTLen)
        return np.dot(self.K, X)

    def mag_cqt(self, x, octaves):
        assert x.size >= self.FFTLen * 2**(octaves - 1)
        result = np.zeros(octaves * self.bins) 
        for o in range(octaves)[::-1]:
            t_oneOct = self.cqt_oneOct(x[-self.FFTLen:])
            for b in range(self.bins):
                result_bin = o*self.bins + b
                result[result_bin] = np.sum(
                    abs(t_oneOct[b*self.atomNr:(b+1)*self.atomNr])
                    ) / self.atomNr
            x = filtfilt(self.B, self.A, x)
            x = x[::2]
        return result

    def best_scale(self, mag_cqt, octaves):
        major = [1,0,1,0,1,1,0,1,0,1,0,1] * octaves
        major_rotations = np.array([major[-p:]+major[:-p] for p in range(12)])
        scale_scores = np.dot(major_rotations, mag_cqt)
        best_scale = np.argmax(scale_scores)
        return major_rotations[best_scale] * (0.25 + 5 * mag_cqt) 


#============================================================================
# Create the Chaco plot.
#============================================================================

def _create_plot_component(obj, cqtkernel):
    # Scale plot
    scale_data = zeros((OCTAVES * BINS, FRAMES_PER_PLOT))
    obj.scale_plotdata = ArrayPlotData()
    obj.scale_plotdata.set_data('scaleimagedata', scale_data)
    scale_plot = Plot(obj.scale_plotdata)
    max_time = float(FRAMES_PER_PLOT * 2**(OCTAVES-1) * cqtkernel.FFTLen) / SAMPLING_RATE
    max_note_index = OCTAVES * BINS
    scale_plot.img_plot('scaleimagedata', 
                        name = 'Scale',
                        xbounds=(0, max_time),
                        ybounds=(0, max_note_index),
                        colormap=hot,
                        )
    scale_range = scale_plot.plots['Scale'][0].value_mapper.range
    scale_range.high = 5
    scale_range.low = 0.0
    scale_plot.title = 'Scale'
    obj.scale_plot = scale_plot

    # CQT plot
    cqt_data = zeros((OCTAVES * BINS, FRAMES_PER_PLOT))
    obj.cqt_plotdata = ArrayPlotData()
    obj.cqt_plotdata.set_data('cqtimagedata', cqt_data)
    cqt_plot = Plot(obj.cqt_plotdata)
    max_time = float(FRAMES_PER_PLOT * 2**(OCTAVES-1) * cqtkernel.FFTLen) / SAMPLING_RATE
    max_note_index = OCTAVES * BINS
    cqt_plot.img_plot('cqtimagedata', 
                        name = 'CQT',
                        xbounds=(0, max_time),
                        ybounds=(0, max_note_index),
                        colormap=hot,
                        )
    cqt_range = cqt_plot.plots['CQT'][0].value_mapper.range
    cqt_range.high = 5
    cqt_range.low = 0.0
    cqt_plot.title = 'CQT'
    obj.cqt_plot = cqt_plot

    container = HPlotContainer()
    container.add(obj.scale_plot)
    container.add(obj.cqt_plot)

    return container

_stream = None

def get_audio_data():
    global _stream
    if _stream is None:
        pa = pyaudio.PyAudio()
        _stream = pa.open(format=pyaudio.paInt16, channels=1, rate=SAMPLING_RATE,
                     input=True, frames_per_buffer=2**(OCTAVES-1) * popup.cqtkernel.FFTLen)
    audio_data  = fromstring(_stream.read(2**(OCTAVES-1) * popup.cqtkernel.FFTLen), dtype=short)
    normalized_data = audio_data / 32768.0
    cqt = popup.cqtkernel.mag_cqt(normalized_data, OCTAVES)
    scale = popup.cqtkernel.best_scale(cqt, OCTAVES)
    return (scale, cqt)

# HasTraits class that supplies the callable for the timer event.
class TimerController(HasTraits):

    def onTimer(self, *args):
        scale, cqt = get_audio_data()
        scale_data = self.scale_plotdata.get_data('scaleimagedata')
        scale_data = hstack((scale_data[:,1:],
                             transpose([scale])))
        self.scale_plotdata.set_data('scaleimagedata', scale_data)

        cqt_data = self.cqt_plotdata.get_data('cqtimagedata')
        cqt_data = hstack((cqt_data[:,1:],
                             transpose([cqt])))
        self.cqt_plotdata.set_data('cqtimagedata', cqt_data)

        return

#============================================================================
# Attributes to use for the plot view.
size = (900,500)
title = "Audio Spectrum"

#============================================================================
# Demo class that is used by the demo.py application.
#============================================================================

class DemoHandler(Handler):

    def closed(self, info, is_ok):
        """ Handles a dialog-based user interface being closed by the user.
        Overridden here to stop the timer once the window is destroyed.
        """

        info.object.timer.Stop()
        return

class Demo(HasTraits):

    plot = Instance(Component)

    controller = Instance(TimerController, ())

    timer = Instance(Timer)

    traits_view = View(
                    Group(
                        Item('plot', editor=ComponentEditor(size=size),
                             show_label=False),
                        orientation = "vertical"),
                    resizable=True, title=title,
                    width=size[0], height=size[1],
                    handler=DemoHandler
                    )

    def __init__(self, **traits):
        super(Demo, self).__init__(**traits)
        self.cqtkernel = CQTKernel(TOP_MIDINUMBER, BINS, SAMPLING_RATE)
        self.plot = _create_plot_component(self.controller, self.cqtkernel)

    def edit_traits(self, *args, **kws):
        self.timer = Timer(20, self.controller.onTimer)
        return super(Demo, self).edit_traits(*args, **kws)

    def configure_traits(self, *args, **kws):
        self.timer = Timer(20, self.controller.onTimer)
        return super(Demo, self).configure_traits(*args, **kws)

popup = Demo()

if __name__ == "__main__":
    try:
        popup.configure_traits()
    finally:
        if _stream is not None:
            _stream.close()

