import nrrd 
import vtk 
import numpy as np 
import numpy.ma as ma
import argparse 
import socket
import logging as log
import os
import sys
import subprocess as sp
from PIL import Image
import datetime

now = datetime.datetime.now()

def tokenize(args):
    single_str = ' '.join(arg for arg in args)
    tokens = single_str.split()
    return tokens

def tile(data, axis):
    if data.shape[axis] < 25:
        print('data is too shallow for tiling along chosen dimension')
        return None
    if axis == 0:
        axes = [0, 1, 2]
        subdata = data[0:25,:,:]
    elif axis == 1:
        axes = [1, 0, 2]
        subdata = data[:,0:25,:]
    else:
        axes = [2, 0, 1]
        subdata = data[:,:,0:25]
    # create a 5 x 5 tiled view of the chosen dimension
    splitted = np.split(subdata, 25, axis=axis)
    splitted = [ s.reshape(data.shape[axes[1]], data.shape[axes[2]]) for s in splitted ]
    rows = [ splitted[5*n:5*(n+1)] for n in range(5) ]
    rows = [ np.concatenate(r, axis=1) for r in rows ]
    square = np.concatenate(rows, axis=0)
    return square

def to_image(data):
    _min = np.min(data)
    _max = np.max(data)
    _range = _max - _min 
    scaled = (data - _min) / _range 
    return Image.fromarray(scaled.astype('uint8'))

def execute(what, cmd):
    cmd = tokenize(cmd)
    log.info(what)
    start_time = datetime.datetime.now()
    log.info(f'started at {start_time}')
    log.info(' '.join(item for item in cmdline))
    completed = sp.run(cmd)
    end_time = datetime.datetime.now()
    log.info(f'ended at {end_time}. duration: {end_time-start_time}\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Determine best ridge extraction scale and visualize resulting measures')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input scalar volume (NRRD)')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output basename')
    parser.add_argument('--log', type=str, default='didthis', help='Basename of log file')
    parser.add_argument('--nctrlpts', type=int, default=5, help='Number of control points for scale space interpolation')
    parser.add_argument('--nscales', type=int, default=11, help='Number of scales to sample')
    parser.add_argument('--maxscale', type=int, default=10, help='Max scale to consider')
    parser.add_argument('-v', '--verbose')
    parser.add_argument('--vis', action='store_true', help='Display samples of computed quantities')
    parser.add_argument('--tiled', type=int, default=0, help='Dimension along which to tile the volumes for display')

    args = parser.parse_args()

    indir = os.path.dirname(args.input)
    outdir = os.path.dirname(args.output)
    base, ext = os.path.splitext(os.path.basename(args.input))
    baseout, extout = os.path.splitext(os.path.basename(args.output))
    baseout = os.path.join(outdir, baseout)
    axes = [0, 0, 0]
    if args.tiled == 0:
        axes = [0, 1, 2]
    elif args.tiled == 1:
        axes = [1, 0, 2]
    else:
        axes = [2, 0, 1]

    # Useful strings
    teem_scsp_base_args = [
        'vprobe',
        f'-i {args.input}',
        '-k scalar -s 1 1 1',
        f'-ssn {args.nctrlpts}',
        f'-ssr 0 {args.maxscale}'
    ]
    teem_scsp_init_args = teem_scsp_base_args + [ f'-sssf blur=%u.nrrd' ]
    teem_scsp_intp_args = teem_scsp_base_args + [ f'-ssrf blur=%u.nrrd', '-ssnd' ]

    logname = args.log + outdir + '.log' 

    log.basicConfig(filename=logname, level=log.DEBUG, format='%(levelname)s:%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    log.info('--------------------------------------------------------')
    log.info('scale_space_navigator.py:')
    log.info('hostname: ' + socket.gethostname())
    log.info('time    : ' + now.strftime('%m/%d/%Y'))
    log.info('--------------------------------------------------------')

    # first, log the current parameter settings
    params = [ f'input: {args.input}', 
               f'output: {args.output}', 
               f'logfile: {args.log}',
               f'nctrlpts: {args.nctrlpts}',
               f'nscales: {args.nscales}',
               f'maxscale: {args.maxscale}',
               f'vis: {args.vis}',
               f'tiled: {args.tiled}' ]

    log.info('Input parameters:')
    for p in params:
        log.info(p)

    # Compute zero-scale values
    cmdline = [ f'vprobe -i {args.input} -k scalar -s 1 1 1 -q gv -o {baseout}_scale=0.0_gv.nrrd' ]
    execute(f'Computing gradient at scale 0', cmdline)
    cmdline = [ f'vprobe -i {args.input} -k scalar -s 1 1 1 -q hess -o {baseout}_scale=0.0_hess.nrrd' ]
    execute(f'Computing Hessian at scale 0', cmdline)
    cmdline = [ f'vprobe -i {args.input} -k scalar -s 1 1 1 -q heval2 -o {baseout}_scale=0.0_heval2.nrrd' ]
    execute(f'Computing ridge strength at scale 0', cmdline)
    cmdline = [ f'cp {args.input} {baseout}_scale=0.0_v.nrrd' ]
    execute('Copying input value', cmdline)

    # Next, compute the control points for Teem scale space interpolation
    cmdline = teem_scsp_init_args + [
        '-q gv',
        f'-o tmp.nrrd' 
    ]
    execute('Computing control points of scale space interpolation', cmdline)

    # Now sample scale space using pre-computed blurred volumes
    lmins = []
    log.info('Computing ridge strength across scale space')
    start_time = datetime.datetime.now()
    log.info(f'started at {start_time}')
    for scale in np.linspace(0, args.maxscale, num=args.nscales, endpoint=True):
        log.info(f'scale: {scale}')
        if scale == 0.0:
            log.info('done')
            continue
        cmdline = teem_scsp_intp_args + [
            '-q heval2',
            f'-ssw {scale}',
            f'-o {baseout}_scale={scale}_heval2.nrrd'
        ]
        execute(f'Ridge strength at scale {scale}', cmdline)
        data, _ = nrrd.read(f'{baseout}_scale={scale}_heval2.nrrd')
        lmins.append(data)
        cmdline = teem_scsp_intp_args + [ 
            '-q gv',
            f'-ssw {scale}',
            f'-o {baseout}_scale={scale}_gv.nrrd'
        ]
        execute(f'Gradient vector at scale {scale}', cmdline)
        cmdline = teem_scsp_intp_args + [ 
            '-q hess',
            f'-ssw {scale}',
            f'-o {baseout}_scale={scale}_hess.nrrd'
        ]
        execute(f'Hessian at scale {scale}', cmdline)
        cmdline = teem_scsp_intp_args + [ 
            '-q v',
            f'-ssw {scale}',
            f'-o {baseout}_scale={scale}_v.nrrd'
        ]
        execute(f'Hessian at scale {scale}', cmdline)

    end_time = datetime.datetime.now()
    log.info(f'ended at {end_time}. duration: {end_time-start_time}\n')

    (w, h, d) = lmins[0].shape

    # Stack ridge strength volumes along scale axis
    log.info('Computing best scale by projection')
    lmins = np.stack(lmins, axis=3)
    bestscales = np.argmin(lmins, axis=3)
    beststr = np.min(lmins, axis=3)    # best lmin across scales
    nrrd.write(f'{baseout}_teem_best_scale.nrrd', bestscales, compression_level=2)
    tiled = tile(bestscales, axis=args.tiled)
    if tiled is not None:
        img = to_image(tiled)
        img.show()

    bestv = np.zeros((w, h, d))     # best values across scales
    bestgv = np.zeros((w, h, d, 3))    # best gradients across scales
    besthess = np.zeros((w, h, d, 9))    # best hessians across scales
    for scale in np.linspace(0, args.maxscale, num=args.nscales, endpoint=True):
        print(f'processing scale {scale}...')
        mask = (bestscales == scale)
        gmask = np.stack([mask, mask, mask], axis=3)
        print(f'gmask.shape: {gmask.shape}')
        hmask = np.stack([mask, mask, mask, mask, mask, mask, mask, mask, mask], axis=3)
        print(f'hmask.shape: {hmask.shape}')
        h, _ = nrrd.read(f'{baseout}_scale={scale}_hess.nrrd')
        print(f'h.shape: {h.shape}')
        h = np.moveaxis(h, 0, -1)
        print(f'now, h.shape: {h.shape}')
        g, _ = nrrd.read(f'{baseout}_scale={scale}_gv.nrrd')
        print(f'g.shape: {g.shape}')
        g = np.moveaxis(g, 0, -1)
        print(f'now, g.shape: {g.shape}')
        v, _ = nrrd.read(f'{baseout}_scale={scale}_v.nrrd')
        hm = ma.array(h, mask=hmask, fill_value=0)
        gm = ma.array(g, mask=gmask, fill_value=0)
        vm = ma.array(v, mask=mask, fill_value=0)
        bestgv = bestgv + gm 
        besthess = besthess + hm
        bestv = bestv + vm

    nrrd.write(f'{baseout}_teem_scsp_value.nrrd', bestv, compression_level=2)
    nrrd.write(f'{baseout}_teem_scsp_gradient.nrrd', bestgv, compression_level=2)
    nrrd.write(f'{baseout}_teem_scsp_hessian.nrrd', besthess, compression_level=2)
    




