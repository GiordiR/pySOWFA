import numpy as np
import scipy.fftpack as scfft
import scipy.signal as scsig


def xcorr_fft(var1, var2=None, maxlags=None, norm='coeff'):
    '''
    cross-correlation using scipy.fftconvolve
    if var1 = var2 compute auto-correlation
    '''
    N = len(var1)
    if var2 is None:
        var2 = var1

    assert len(var1) == len(var2)

    if maxlags is None:
        maxlags = N - 1
        lags = np.arange(0, 2*N-1)
    else:
        assert maxlags < N
        lags = np.arange(N-maxlags-1, N+maxlags)

    res = scsig.fftconvolve(var1, var2[::-1], mode='full')

    Nf = float(N)
    if norm == 'biased':
        res = res[lags] / Nf
    elif norm == 'coeff':
        var1_rms = np.sqrt(np.mean(var1**2))
        var2_rms = np.sqrt(np.mean(var2 ** 2))
        rms = var1_rms * var2_rms
        res = res[lags] / rms / Nf
        res = res[round((len(res)-1)/2):-1]

    lags = np.arange(0, maxlags)

    return res, lags

def FFT(sig, samplefrq, nparseg=None, detrend='constant'):
    '''
    Estimate power spectral density using Welch's method
    '''

    # old crap...
    # siglength = sig.shape[0]
    # NFFT = nextpow2(siglength)
    ##    NFFT = nextpow2(sig.shape(-1,)[0])
    # amp = spfft.fft(sig.reshape(-1,),NFFT)/siglength
    # frq = samplefrq/2*np.linspace(0,1,NFFT/2+1)
    ##    frq = spfft.fftfreq(siglength, (samplefrq)**-1)
    ##    frq = frq[:NFFT/2+1]
    # amp = 2*abs(amp[:NFFT/2+1])
    # return frq,amp

    frq, psd = scsig.welch(sig, fs=samplefrq, window='hann',
                           noverlap=None, nperseg=nparseg,
                           return_onesided=True, scaling='density',
                           axis=-1, detrend=detrend)
    return frq, psd
