from __future__ import print_function
from astropy.convolution import convolve_fft, Gaussian2DKernel
from radio_beam import Beam
from astropy.io import fits
import astropy.units as u
import scipy
import sys
from auxcodes import report

def do_convolve(filename,resolution,outfile=None):
    report('Reading files')
    hdu=fits.open(filename,memmap=False)
    old_b=Beam.from_fits_header(hdu[0].header)
    new_b=Beam(resolution*u.arcsec)
    rr=new_b.sr/old_b.sr
    print('Beam area ratio is %.2f' % rr)
    dcb=new_b.deconvolve(old_b)
    gauss_kern = dcb.as_kernel(1.5*u.arcsec)
    image=hdu[0].data[0,0]
    report('FFT')
    fft_mp = lambda a: scipy.fft.fftn(a, workers=-1)
    ifft_mp = lambda a: scipy.fft.ifftn(a, workers=-1)

    smoothed_data_gauss = convolve_fft(image, gauss_kern, allow_huge=True, fftn=fft_mp, ifftn=ifft_mp)
    report('Scaling')
    hdu[0].data[0,0]=smoothed_data_gauss*rr
    hdu[0].header.update(new_b.to_header_keywords())
    report('Write to disk')
    if outfile is None:
        outfile=filename.replace('.fits','_convolved.fits')
    hdu.writeto(outfile,overwrite=True)

if __name__=='__main__':
    filename=sys.argv[1]
    resolution=float(sys.argv[2])
    do_convolve(filename,resolution)
