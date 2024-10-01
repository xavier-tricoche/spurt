import numpy as np 
import argparse as ap 
import scipy as sp 
import os 
import nrrd 
from PIL import Image, ImageOps 
from matplotlib import pyplot as plt 
import vtk 
from vtk.util import numpy_support as nps
from collections import deque

def fix_header(data, header=None):
    n = data.ndim
    if header is None:
        spacings = np.ones(n)
        sizes = np.array(data.shape)
        return { 'sizes': sizes, 'spacings': spacings }
    else:
        if 'spacings' not in header.keys():
            header['spacings'] = np.ones(n)
        if 'sizes' not in header.keys():
            header['sizes'] = np.array(data.shape)
        return header

def read(filename):
    ext = os.path.splitext(filename)[1].lower()
    if ext == '.nrrd' or ext == '.nhdr':
        d, h = nrrd.read(filename)
        return d, fix_header(d, h)
    elif ext in [ '.png', '.jpg', '.jpeg', '.jp2', '.bmp', '.ppm', '.tif', '.tiff', '.gif']:
        im = np.array(Image.open(filename))
        return im, fix_header(im)
    else:
        raise ValueError(f'Unrecognized image file extension: {ext}')

def export(data, filename):
    basename, ext = os.path.splitext(filename)
    if isinstance(data, vtk.vtkImageData):
        writer = vtk.vtkXMLImageDataWriter()
        filename = basename + '.vti'
    elif isinstance(data, vtk.vtkPolyData):
        writer = vtk.vtkPolyDataWriter()
        filename = basename + '.vtp'
    elif isinstance(data, np.ndarray):
        filename = basename + '.nrrd'
        nrrd.write(filename, data)
        return 
    elif isinstance(data, vtk.vtkDataSet):
        filename = basename + '.vtk'
        writer = vtk.vtkDataSetWriter()
    else:
        raise ValueError(f'Unrecognized dataset type for {data}')
    
    writer.SetFileName(filename)
    writer.SetInputData(data)
    writer.Write()

def gradient(data, header):
    gradx, grady = np.gradient(data, *header['spacings'])
    newh = header.copy()
    newh['sizes'] = np.concatenate((header['sizes'], [2]))
    newh['spacings'] = np.concatenate((header['spacings'], [np.nan]))
    return np.stack([gradx, grady], axis=-1), newh

def hessian(gradient, header):
    width, height = gradient.shape[:2]
    gradx = gradient[:,:,0]
    grady = gradient[:,:,1]
    Hxx, Hxy = np.gradient(gradx, *header['spacings'][:2])
    Hyy = np.gradient(grady, axis=1, *[ header['spacings'][1] ])
    H = np.stack([Hxx, Hxy, Hxy, Hyy], axis=-1)
    H = H.reshape((width, height, 2, 2))
    newh = header.copy()
    sz = header['sizes']
    spc = header['spacings']
    newh['sizes'] = np.array(H.shape)
    newh['spacings'] = np.concatenate((spc, [np.nan]))
    return H, newh

def evals(M):
    return np.linalg.eigvalsh(M)

def compute_derivatives(image, header):
    grad, grad_header = gradient(image, header)
    print('computed gradient')
    hess, hess_header = hessian(grad, grad_header)
    print('computed hessian')
     # min eigenvalue of Hessian is ridge strength
    ridge_strength = evals(hess)[:,:,0]
    print('compute ridge strength')
    return grad, grad_header, hess, hess_header, ridge_strength

def isoline(data, value=0):
    c = vtk.vtkContourFilter()
    c.SetValue(0, value)
    c.SetInputData(data)
    c.Update()
    return c.GetOutput()

def which(apair, done):
    if done[apair[0]]:
        return apair[1]
    else:
        return apair[0]

def neighbor(edge, id0, throw=False):
    if id0 == edge[0]:
        return edge[1]
    elif id0 == edge[1]:
        return edge[0]
    else:
        raise ValueError(f'given id ({id0}) is not part of edge {edge}')

def extend_cc(acc, p2c, segments, done_pts):
    seed = acc[-1]
    cell_ids = p2c[seed]
    if not cell_ids: 
        print(f'Warning: extend_cc given pt {seed} that has no neightbors')
        return acc 
    elif len(cell_ids) == 1:
        neigh = neighbor(segments[cell_ids[0]], seed)
        if not done_pts[neigh]:
            acc.append(neigh)
            done_pts[neigh] = True
        print(f'Warning: extend_cc give pt {seed} that has only one neighbor')
        return acc
    else:
        n1 = neighbor(segments[cell_ids[0]], seed)
        n2 = neighbor(segments[cell_ids[1]], seed)
        if done_pts[n1] and not done_pts[n2]:
            acc.append(n2)
            done_pts[n2] = True 
            return extend_cc(acc, p2c, segments, done_pts)
        elif not done_pts[n1] and done_pts[n2]:
            acc.append(n1)
            done_pts[n1] = True 
            return extend_cc(acc, p2c, segments, done_pts)
        elif not done_pts[n1] and not done_pts[n2]:
            raise ValueError(f'forward extension from {seed} reached two unprocessed vertices {n1} and {n2}')
        else:
            print(f'warning: extend_cc at {seed} cannot proceed')
            return acc
        
def get_cell_ids(cells):
    ids_ = vtk.reference([])
    n_ = vtk.reference(0)
    cells.GetNextCell(n_, ids_)
    n = n_.get()
    ids = ids_.get()
    return ids

def add_or_insert(key, value, adict):
    if key in adict.keys():
        adict[key].append(value)
    else:
        adict[key] = [value]

def cc(contour, minlength=3):
    segments = contour.GetLines()
    npts = contour.GetNumberOfPoints()
    print(f'there are {segments.GetNumberOfCells()} in cc input')
    p2e = {}
    cells = []
    # compute point to edge(s) mapping
    segments.InitTraversal()
    for i in range(segments.GetNumberOfCells()):
        ids = get_cell_ids(segments)
        cells.append(ids)
        for id in ids:
            add_or_insert(id, i, p2e)
    ncells = len(cells)

    # now compute connected components
    done_pts = np.zeros((npts), dtype=bool)
    newcells = []
    for ptid, cellids in p2e.items():
        print(f'processing vertex #{ptid}')
        # skip current point if already seen
        if done_pts[ptid]: 
            print('already seen: skip')
            continue 
        done_pts[ptid] = True
        # skip also if it belongs to 0 or 1 segment
        if not cellids: 
            print(f'no cell: skip')
            continue
        elif len(cellids) == 1:
            id = neighbor(cells[cellids[0]], ptid)
            print(f'neighbor({ptid}) = {id}')
            if done_pts[id]:
                print(f'Warning: neighbor of current seed has already been processed')
                continue
            else:
                done_pts[id] = True
                fwd = extend_cc([id], p2e, cells, done_pts)
                if len(fwd) >= minlength:
                    newcells.append(fwd)
        else:
            id1, id2 = cellids
            leftid = neighbor(cells[id1], ptid)
            rightid = neighbor(cells[id2], ptid)
            if done_pts[leftid] or done_pts[rightid]:
                raise ValueError(f'ERROR in CC: point {ptid} is connected to  vertices already seen (either {leftid} or {rightid})')
            acell = deque([leftid, ptid, rightid])
            done_pts[leftid] = True 
            done_pts[rightid] = True
            fwd = extend_cc([rightid], p2e, cells, done_pts)
            bwd = extend_cc([leftid], p2e, cells, done_pts)
            acell.extend(fwd)
            acell.extendleft(bwd)
            if len(acell) >= minlength:
                newcells.append(list(acell))

    cell_array = vtk.vtkCellArray()
    for c in newcells:
        cell_array.InsertNextCell(len(c))
        for p in c:
            cell_array.InsertCellPoint(p)

    contour.SetLines(cell_array)
    return contour

def is_ridge_by_curvature(H, g, Hg):
    if np.linalg.norm(g) == 0.:
        # compute determinant and trace of H
        if np.det(H) > 0:
            return True
    else:
        c = np.array([-g[1], g[0]])
        Hc = np.matmul(H, c)
        if abs(g[0]) > abs(g[1]):
            lg = Hg[0]/g[0]
            lc = Hc[1]/c[1]
        else:
            lg = Hg[1]/g[1]
            lc = Hc[0]/c[0]
        return lc < 0

def remove_valleys(dataset, minsize=2):
    npts = dataset.GetNumberOfPoints()
    is_ridge = np.zeros((npts), dtype=bool)
    for i in range(npts):
        p = dataset.GetPoint(i)
        g = dataset.GetPointData().GetArray('g').GetTuple(i)
        Hg = dataset.GetPointData().GetArray('Hg').GetTuple(i)
        flatH = dataset.GetPointData().GetArray('H').GetTuple(i)
        H = np.array([[flatH[0], flatH[1]], [flatH[1], flatH[2]]])
        is_ridge[i] = is_ridge_by_curvature(H, g, Hg)
    
    nptsremoved = npts - np.sum(is_ridge)
    print(f'{nptsremoved} crease points were removed ({float(nptsremoved)/float(npts)*100.}%)')
    old_lines = dataset.GetLines()
    nlines = old_lines.GetNumberOfCells()

    print(f'there are {nlines} old lines to go through')
    
    old_lines.InitTraversal()
    lines_array = []
    stopped = True
    for c in range(nlines):
        ids = get_cell_ids(old_lines)
        for id in ids:
            if not is_ridge[id]:
                stopped = True
            elif stopped:
                lines_array.append([id])
                stopped = False 
            else:
                lines_array[-1].append(id)
    new_lines = vtk.vtkCellArray()
    lengths = []
    for l in lines_array:
        if len(l) >= minsize:
            new_lines.InsertNextCell(len(l))
            for id in l:
                new_lines.InsertCellPoint(id)
            lengths.append(len(l))

    lengths = np.array(lengths)
    print(f'after filtering cells based on ridge points there are {new_lines.GetNumberOfCells()} lines left')
    print(f'Lengths of cells range from {np.min(lengths)} to {np.max(lengths)}')
    print(f'There are {np.sum(lengths >= 10)} lines longer than 10')

    '''
    connect = vtk.vtkConnectivityFilter()
    connect.SetInputData(dataset)
    connect.SetExtractionModeToAllRegions()
    connect.ColorRegionsOn()
    connect.Update()
    ccs = connect.GetOutput()
    print('connected components:\n', ccs)
    regionids = ccs.GetPointData().GetArray('RegionId')
    regionids = nps.vtk_to_numpy(regionids)
    maxid = np.max(regionids)
    nregions = np.unique(regionids).shape[0]
    print(f'Max region id is {maxid}, number of regions is {nregions}')
    sortids = np.argsort(regionids)
    regions = {}
    for i, r in enumerate(regionids):
        if r in regions.keys():
            regions[r].append(i)
        else:
            regions[r] = [i]
    sizes = []
    for r, ids in regions.items():
        sizes.append(len(ids))
    _min = np.min(sizes)
    _max = np.max(sizes)
    _mean = np.mean(sizes)
    sizes = np.sort(sizes)
    print(f'cc sizes: min: {_min}, mean: {_mean}, max: {_max}, largest regions: {sizes[-10:-1]}')
    '''
    
    new_dataset = vtk.vtkPolyData()
    newpts = vtk.vtkPoints()
    newpts.DeepCopy(dataset.GetPoints())
    new_dataset.SetPoints(newpts)
    new_dataset.SetLines(new_lines)
    new_dataset.GetPointData().DeepCopy(dataset.GetPointData())
    # return cc(new_dataset)
    return new_dataset

def marching_ridges(_image, _header, threshold, minlength, upsampling=1):
    from image import DifferentiableImage as DImage
    dataset = DImage(_image, _header)
    if upsampling > 1:
        nx = DImage.



    width, height = header['sizes'][:2]
    grad, grad_header, \
    hess, hess_header, \
    ridge_strength = compute_derivatives(image, header)
    evals, evecs = np.linalg.eigh(hess)
    lmin = evals[:,:,0] 
    emin = evecs[:,:,:,0]
    lines = {}
    for i in range(width-1):
        for j in range(height-1):
            ids = [ [i, j], [i+1, j], [i+1, j+1], [i, j+1] ]
            ls = [ lmin[a, b] for a, b in ids ]
            if np.min(ls) > threshold: continue

            emins = [ emin[a, b, :] for a, b in ids ]
            C = np.zeros((2,2), dype=float)
            for e in emins:
                C += np.outer(e, e)
            C /= 4.
            _, evs = np.eigh(C)
            eavg = evs[:,0]
            for e in emins:
                dot = np.dot(e, eavg)
                if dot < 0:
                    e *= -1
            for edge in range(4):
                id0 = ids[edge]
                id1 = ids[(edge+1)%4]




def peikert_sadlo(image, header, threshold, minlength=3):
    width, height = header['sizes'][:2]
    grad, grad_header, \
    hess, hess_header, \
    ridge_strength = compute_derivatives(image, header)
    hess_grad = np.einsum('...ij,...j', hess, grad)
    export(hess, 'hessian.nrrd')
    export(grad, 'gradient.nrrd')
    export(hess_grad, 'hessian_times_gradient.nrrd')
    export(ridge_strength, 'ridge_strength.nrrd')

    ps_criterion = np.cross(hess_grad, grad) # 
    print('computed criterion value')
    export(ps_criterion, 'ps_criterion.nrrd')

    ps_criterion_dataset = vtk.vtkImageData()
    ps_criterion_dataset.SetDimensions(width, height, 1)
    ps_criterion_dataset.GetPointData().SetScalars(nps.numpy_to_vtk(ps_criterion.ravel(order='F')))
    export(ps_criterion_dataset, 'ps_criterion.vti')

    # superset of ridge lines as zero level set of the criterio
    raw = isoline(ps_criterion_dataset, 0) 
    histo, bins = np.histogram(ridge_strength, bins=20)
    txt = ''
    for i, val in enumerate(histo):
        txt += f'({bins[i]}) {val} '
    txt += f'{bins[-1]}'
    print(txt)
    print(f'Strength percentiles: 10%: {np.percentile(ridge_strength, 10)}, 5%: {np.percentile(ridge_strength, 5)}, 2%: {np.percentile(ridge_strength, 2)}, 1%: {np.percentile(ridge_strength, 1)}, 0.5%: {np.percentile(ridge_strength, 0.5)}')

    image = vtk.vtkImageData()
    image.SetDimensions(width, height, 1)
    vtk_ridge_strengtharray = nps.numpy_to_vtk(ridge_strength.ravel(order='F'))
    vtk_ridge_strengtharray.SetName('ridge strength')
    image.GetPointData().SetScalars(vtk_ridge_strengtharray)

    fake_data = np.zeros((width, height))
    vals = np.linspace(0, height, height)
    for h in range(height):
        fake_data[:,h] = np.linspace(vals[h], vals[h]+1, width)
    print(fake_data)
    vtk_fake_data = nps.numpy_to_vtk(fake_data.reshape((width*height), order='F'))
    fake_image = vtk.vtkImageData()
    fake_image.SetDimensions(width, height, 1)
    fake_image.GetPointData().SetScalars(vtk_fake_data)
    export(fake_image, 'fake.vti')

    flattened_hess = np.delete(hess.reshape((width*height, 4), order='F'), 2, axis=1)
    vtk_hess_array = nps.numpy_to_vtk(flattened_hess)
    vtk_hess_array.SetName('H')
    export(flattened_hess, 'flattened_hessian.vti')
    image.GetPointData().AddArray(vtk_hess_array)

    vtk_hess_grad_array = nps.numpy_to_vtk(hess_grad.reshape((width*height, 2), order='F'))
    vtk_hess_grad_array.SetName('Hg')
    image.GetPointData().AddArray(vtk_hess_grad_array)

    vtk_grad_array = nps.numpy_to_vtk(grad.reshape((width*height, 2), order='F'))
    vtk_grad_array.SetName('g')
    image.GetPointData().AddArray(vtk_grad_array)

    probe = vtk.vtkProbeFilter()
    probe.SetInputData(raw)
    probe.SetSourceData(image)
    probe.Update()
    raw_plus_values = probe.GetOutput()
    print('Raw\n', raw_plus_values)

    clip = vtk.vtkClipPolyData()
    clip.InsideOutOn()
    clip.SetValue(threshold)
    clip.SetInputData(raw_plus_values)
    clip.Update()
    crease_lines = clip.GetOutput()
    print('Crease lines\n', crease_lines)

    ridge_lines = remove_valleys(crease_lines, minlength)
    print(f'after cc computation, there are {ridge_lines.GetLines().GetNumberOfCells()} lines')

    export(image, 'image.vti')
    export(probe.GetOutput(), 'isolines.vtp')
    export(crease_lines, 'crease_lines.vtp')
    export(ridge_lines, 'ridge_lines.vtp')

if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Compute ridge lines in 2D images using Peikert and Sadlo\'s method presented in their 2008 paper.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input image')
    parser.add_argument('-t', '--threshold', type=float, default=0., help='Threshold on ridge strength (negative) in raw ridge filtering')
    parser.add_argument('--valley', action='store_true', help='Extract valley lines instead of ridges')
    parser.add_argument('--blur', type=float, default=0, help='Perform Gaussian blur prior to extracting ridges')
    parser.add_argument('--length', type=int, default=2, help='Minimum required ridge line length')

    args = parser.parse_args()
    image, header = read(args.input)
    if args.valley:
        image *= -1.
    if args.blur > 0:
        image = sp.ndimage.gaussian_filter(image, sigma=args.blur)

    peikert_sadlo(image=image, header=header, threshold=args.threshold, minlength=args.length)