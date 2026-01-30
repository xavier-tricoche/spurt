import scipy as sp 
import nrrd 
import numpy as np
import argparse as ap 
import time
import os
from scipy.ndimage import gaussian_filter
from tqdm import tqdm


class Clock:
    def __init__(self):
        self.t = time.time()
    
    def tick(self):
        self.t = time.time()
    
    def tock(self):
        return time.time() - self.t

def value(image, sigma):
    return gaussian_filter(image, sigma=sigma, order=0, mode='nearest')

def gradient(image, sigma):
    w, h, d = image.shape
    gv0 = gaussian_filter(image, sigma=sigma, order=(1,0,0), mode='nearest')
    gv1 = gaussian_filter(image, sigma=sigma, order=(0,1,0), mode='nearest')
    gv2 = gaussian_filter(image, sigma=sigma, order=(0,0,1), mode='nearest')
    _gradient = np.zeros((w, h, d, 3))
    _gradient[...,0] = gv0
    _gradient[...,1] = gv1
    _gradient[...,2] = gv2
    _gradient *= sigma
    return _gradient

def hessian(image, sigma):
    w, h, d = image.shape
    h00 = gaussian_filter(image, sigma=sigma, order=(2,0,0), mode='nearest')
    h01 = gaussian_filter(image, sigma=sigma, order=(1,1,0), mode='nearest')
    h02 = gaussian_filter(image, sigma=sigma, order=(1,0,1), mode='nearest')
    h11 = gaussian_filter(image, sigma=sigma, order=(0,2,0), mode='nearest')
    h12 = gaussian_filter(image, sigma=sigma, order=(0,1,1), mode='nearest')
    h22 = gaussian_filter(image, sigma=sigma, order=(0,0,2), mode='nearest')
    _hessian = np.zeros((w, h, d, 3, 3))
    _hessian[...,0,0] = h00
    _hessian[...,0,1] = h01
    _hessian[...,0,2] = h02
    _hessian[...,1,0] = h01
    _hessian[...,1,1] = h11
    _hessian[...,1,2] = h12
    _hessian[...,2,0] = h02
    _hessian[...,2,1] = h12
    _hessian[...,2,2] = h22
    _hessian *= sigma*sigma
    return _hessian

def ridge_strength(hessian):
    return np.linalg.eigvalsh(hessian)[...,0]
    
def sample_scale(img, sigma):
    print(f'Compute image derivatives: sigma={sigma}')
    m, n, p = img.shape
    clock = Clock()
    v = value(img, sigma)
    print(f'blurred computed in {clock.tock()} seconds')
    clock.tick()
    gv = gradient(img, sigma)
    print(f'Gradient computed in {clock.tock()} seconds')
    clock.tick()
    hess = hessian(img, sigma)
    print(f'Hessian computed in {clock.tock()} seconds')
    clock.tick()
    heval2 = ridge_strength(hess)
    print(f'Ridge strength computed in {clock.tock()} seconds')
    return v, gv, hess, heval2

if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Extract/highlight ridges from 3D images using scikit-image')
    parser.add_argument('-i', '--input', type=str, required=True, help='Name of input image file')
    parser.add_argument('-o', '--output', type=str, help='Basename of output image files')
    parser.add_argument('-s', '--sigma', type=float, nargs=3, default=[1,8,0.5], help='Sigma coefficients for Gaussian blurring')
    args = parser.parse_args()

    scales = np.arange(args.sigma[0], args.sigma[1], step=args.sigma[2])
    data, _ = nrrd.read(args.input)
    w, h, d = data.shape
    best_value = np.zeros((w, h, d))
    best_hessian = np.zeros((w, h, d, 3, 3))
    best_emin = np.ones((w, h, d))
    best_scale = np.zeros((w, h, d))
    best_gradient = np.zeros((w, h, d, 3))
    for s in tqdm(scales):
        start = time.time()
        v, gv, hess, emin = sample_scale(data, s)
        new_best = emin<best_emin
        best_scale[new_best] = s
        best_hessian[new_best,:,:] = hess[new_best,:,:]
        best_emin[new_best] = emin[new_best]
        best_gradient[new_best,:] = gv[new_best,:]
        best_value[new_best] = v[new_best]
        end = time.time()
        print(f'Scale processing at sigma={s} took {end-start} s.')

    basename = os.path.splitext(args.output)[0]
    print(best_hessian.shape)
    print(best_gradient.shape)

    nrrd.write(basename + '_scalespace_value.nrrd', best_value, compression_level=3)
    nrrd.write(basename + '_scalespace_scale.nrrd', best_scale, compression_level=3)
    nrrd.write(basename + '_scalespace_hessian.nrrd', best_hessian, compression_level=3)
    nrrd.write(basename + '_scalespace_strength.nrrd', best_emin, compression_level=3)
    nrrd.write(basename + '_scalespace_gradient.nrrd', best_gradient, compression_level=3)

    