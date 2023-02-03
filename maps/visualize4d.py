import vtk
import numpy as np
import argparse
import nrrd
import sys
import os
import pandas as pd
from scipy import spatial
import random
from tqdm import tqdm
import re
import colors
from sklearn import manifold
import clipping_plane as clip

from vtk.util.numpy_support import *

with_axes = False

name_to_column = {
    'p1': 'xd',
    'p2': 'yd',
    'q1': 'x',
    'q2': 'y',
    'r' : 'norm'
}

def draw_axes(renderer, inter, dims, scale=0.25, offset=(0,0,0)):
    print(f'scale={scale}')
    om = vtk.vtkOrientationMarkerWidget()
    cube = colors.make_cube_axis_actor(dims)
    om.SetOrientationMarker(cube)
    om.SetInteractor(inter)
    om.SetDefaultRenderer(renderer)
    return om

def mindist(points, nneighs=2):
    kdtree = spatial.KDTree(points)
    dist, ids = kdtree.query(points, k=nneighs+1)
    dist = dist[:,1:]
    _max = np.max(dist)
    _mean = np.mean(dist)
    # print('max 2nd min distance is {}'.format(max))
    return [_mean, _max]

def make_glyphs(dataset, scale=False, scale_factor=0.001, sres=6):
    if sres >= 5:
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(scale_factor)
        sphere.SetThetaResolution(sres)
        sphere.SetPhiResolution(sres)
        glyphs = vtk.vtkGlyph3D()
        glyphs.SetSourceConnection(sphere.GetOutputPort())
        glyphs.SetInputData(dataset)
        glyphs.ScalingOff()
        glyphs.SetColorModeToColorByScalar()
        return glyphs
    else:
        n = dataset.GetNumberOfPoints()
        ids = vtk.vtkCellArray()
        ids.InitTraversal()
        for i in range(n):
            ids.InsertNextCell(1)
            ids.InsertCellPoint(i)
        dataset.SetVerts(ids)
        alg = vtk.vtkTrivialProducer()
        alg.SetOutput(dataset)
        return alg

def find_files_(expression):
    output = []
    path, template = os.path.split(expression)
    path = os.path.expanduser(path)
    print('path={}'.format(path))
    print('template={}'.format(template))
    with os.scandir(path) as filenames:
        # print(f'filenames={filenames}')
        for fname in tqdm(filenames):
            fname = fname.name
            # print(f'fname={fname}')
            m = re.match(template, fname)
            if m is not None:
                output.append(os.path.join(path, fname))
    return output

def find_files(input):
    output = []
    if isinstance(input, list):
        for expression in input:
            output.extend(find_files_(expression))
    else:
        output.extend(find_files_(input))
    # print('result:\n{}'.format(output))
    return output

def fit_surface(df, dims, args):
    ids = df['index'].unique()
    isomap = manifold.Isomap(n_components=dims)
    cells = vtk.vtkCellArray()
    cellindex = []
    for id in ids:
        subdf = df[df['index']==id]
        # print('subdf={}'.format(subdf))
        indices = subdf.index
        # print('row ids are {}'.format(subdf.index))
        flat = isomap.fit_transform(subdf[['x', 'y', 'xd', 'yd']].to_numpy())
        tris = spatial.Delaunay(flat)
        for t in tris.simplices:
            cellindex.append(id)
            cells.InsertNextCell(3)
            for j in range(3):
                cells.InsertCellPoint(indices[t[j]])
    cellindex = numpy_to_vtk(np.array(cellindex))
    cellindex.SetName('index')
    return cells, cellindex

def plot_surfaces(df, dims, args):
    cells, values = fit_surface(df, dims, args)

def bounds(df, dims):
    mins = df[[dims[0], dims[1], dims[2]]].min(axis=1)
    maxs = df[[dims[0], dims[1], dims[2]]].max(axis=1)
    return mins, maxs

def missing_axis(dims):
    for axis in ['x', 'y', 'xd', 'yd']:
        if axis not in dims:
            return axis
    raise ValueError('No missing axis in {}'.format(dims))

def stereographic_projection(df, dims, lw):
    t = missing_axis(dims)
    x = df[dims[0]]/(lw-df[t])
    y = df[dims[1]]/(lw-df[t])
    z = df[dims[2]]/(lw-df[t])
    return numpy_to_vtk(np.array([x, y, z]).transpose())

def create_clouds(df, all_dims, args):
    clouds = []
    # all_dims = [ ['x', 'y', 'xd'], ['x', 'y', 'yd'],
    # ['x', 'xd', 'yd'], ['y', 'xd', 'yd'] ]

    for f, dims in enumerate(all_dims):
        # coords = all_dims[frame]
        if args['flatten'] is None:
            if args['stereo'] is None:
                c = numpy_to_vtk(df[[dims[0], dims[1], dims[2]]].to_numpy())
            else:
                c = stereographic_projection(df, dims, args['stereo'])
                print('c=\n{}'.format(c))
        else:
            axis = missing_axis(dims)
            eps = args['flatten']
            subdf = df[ (df[axis]<eps) & (df[axis]>-eps) ]
            c = numpy_to_vtk(subdf[[dims[0], dims[1], dims[2]]])
        points = vtk.vtkPoints()
        points.SetData(c)
        cloud = vtk.vtkPolyData()
        cloud.SetPoints(points)
        clouds.append(cloud)

    if args['flatten'] is None:
        return clouds, df
    else:
        return clouds, subdf

def import_orbits(input, args):
    names = []
    if isinstance(input, list):
        names = input
    else:
        names.append(input)

    if len(names) == 0:
        print('ERROR: no filenames in input')
        print(f'input={input}')
        sys.exit(0)

    main_df = None
    random.seed(130875)
    if args['total'] is not None:
        if args['max'] is not None and len(names)*args['max']>args['total']:
            nbfiles = int(args['total']/args['max'])
            n = len(names)
            random.shuffle(names)
            names = names[:nbfiles]
            print('we are choosing at random {} out of {}'.format(nbfiles, n))
    for i, name in enumerate(tqdm(names)):
        # print('processing {}'.format(name))
        data, header = nrrd.read(name)
        df = pd.DataFrame(data, index=['x', 'y', 'xd', 'yd']).transpose()
        df['order'] = np.linspace(0, 1, data.shape[1])
        df['index'] = i

        if args['max'] is not None and data.shape[1]>args['max']:
            # print('reducing orbit from {} to {} points'.format(data.shape[1], args['max']))
            df = df.sample(n=args['max'], replace=False, ignore_index=False)

        if main_df is not None:
            main_df = pd.concat([main_df, df], ignore_index=True)
        else:
            main_df = df
    main_df.reset_index()
    return main_df

class Frame:
    def __init__(self):
        self.renderer = vtk.vtkRenderer()
        self.cloud = None
        self.clipper = clip.ClippingPlane()
        self.active_clipper = False
        self.dims = [ 'x', 'y', 'xd' ]
        self.plane_info = { 'origin': [0,0,0], 'normal': [0,0,1] }

    def toggle_clipper(self):
        if self.active_clipper:
            self.clipper.Deactivate()
        else:
            self.clipper.Activate()
        self.active_clipper = not self.active_clipper

class FourDimensionalVisualizer:
    def __init__(self, input, args):
        # something
        self.input = input
        self.args = args
        self.data = import_orbits(input, args)
        self.window = vtk.vtkRenderWindow()
        self.multiren = True
        self.frames = []
        self.interactor = vtk.vtkRenderWindowInteractor()
        self.interactor.SetRenderWindow(self.window)
        self.xmins = [ 0, 0.5, 0, 0.5 ]
        self.ymins = [ 0, 0, 0.5, 0.5 ]
        self.oms = []
        if args['coordinates'] is not None:
            self.multiren = False
            self.all_dims = [ self.args['coordinates'] ]
            frame = Frame()
            frame.dims = self.args['coordinates']
            self.frames.append(frame)
        else:
            self.all_dims = [ ['x', 'y', 'xd'], ['x', 'y', 'yd'],
            ['x', 'xd', 'yd'], ['y', 'xd', 'yd'] ]
            for i in range(4):
                frame = Frame()
                frame.dims = self.all_dims[i]
                frame.renderer.SetViewport(self.xmins[i], self.ymins[i], self.xmins[i]+0.5, self.ymins[i]+0.5)
                self.frames.append(frame)
        self.clippings = []
        self.active_clippings = [ False, False, False, False ]

    def clipper_callback(self):
        for f in self.frames:
            f.toggle_clipper()
        self.window.Render()

    def flatten_callback(self, eps=1.0e-3):
        if self.args['flatten'] is None:
            self.args['flatten'] = eps
        else:
            self.args['flatten'] = None
        self.run(False)

    def stereo_callback(self, increase):
        if self.args['stereo'] is None:
            self.args['stereo'] = 5
        elif increase:
            self.args['stereo'] = self.args['stereo'] + 0.5
        else:
            self.args['stereo'] = self.args['stereo'] - 0.5
        self.run(False)

    def key_pressed_callback(self, obj, event):
        new_key = obj.GetKeySym()
        if new_key.isdigit():
            key = int(new_key)
            print(f'you pressed key {key}')
        elif new_key == 'c':
            print('you requested clipping')
            self.clipper_callback()
        elif new_key == 'f':
            print('you requested flattening')
            self.flatten_callback()
        elif new_key == 'plus':
            print('you requested an increase in light distance')
            self.stereo_callback(True)
        elif new_key == 'minus':
            print('you requested a decrease in light distance')
            self.stereo_callback(False)
        else:
            print('you typed {}'.format(new_key))

    def run(self, first_time=True):
        if first_time:
            if self.args['color'] == 'r' or self.args['color'] == 'norm':
                self.data['norm'] = self.data.apply(np.linalg.norm, axis=1)

        self.clouds, self.visible = create_clouds(self.data, self.all_dims, self.args)
        for i, f in enumerate(self.frames):
            f.cloud = self.clouds[i]

        att_name = None
        do_random = False
        if self.args['color'] in list(self.data.columns):
            att_name = self.args['color']
        elif self.args['color'] in name_to_column.keys():
            att_name = name_to_column[self.args['color']]
        elif self.args['color'] == 'random':
            do_random = True
        else:
            raise ValueError('unrecognized attribute name for color: {}'.format(args['color']))
        if not do_random:
            self.attribute = numpy_to_vtk(self.visible[att_name].to_numpy())
            self.attribute.SetName(att_name)
        else:
            self.attribute = colors.create_vtk_colors(self.visible['index'].to_numpy())
            self.attribute.SetName('random')
        self.cmap = colors.make_colormap('viridis', self.attribute.GetRange())

        if self.args['triangulate']:
            cells, values = fit_surface(self.visible, 3, self.args)

        print(f'the range of values in selected attribute is {self.attribute.GetRange()}')
        print(f'attribute: {self.attribute}')

        print('there are {} frames'.format(len(self.frames)))
        for i, f in enumerate(self.frames):
            if not first_time:
                f.cloud.GetPointData().RemoveArray(att_name)
            f.cloud.GetPointData().AddArray(self.attribute)
            f.cloud.GetPointData().SetActiveScalars(self.attribute.GetName())
            if self.args['triangulate']:
                f.cloud.SetPolys(cells)
                f.cloud.GetCellData().AddArray(values)
                f.cloud.GetCellData().SetActiveScalars('index')
                self.cmap2 = colors.make_colormap('viridis', values.GetRange())

            glyphs = make_glyphs(f.cloud, self.args['scale'], self.args['factor'], self.args['resolution'])

            if first_time:
                f.clipper.Initialize(glyphs, f.renderer, self.interactor)
                f.clipper.mapper.ScalarVisibilityOn()
                f.clipper.actor.GetProperty().RenderPointsAsSpheresOn()
                f.clipper.actor.GetProperty().SetPointSize(10)
            else:
                f.clipper.Update(glyphs)
            f.clipper.mapper.SetLookupTable(self.cmap)

            if self.args['triangulate']:
                m2 = vtk.vtkPolyDataMapper()
                m2.SetInputData(f.cloud)
                m2.ScalarVisibilityOn()
                m2.SetLookupTable(cmap2)
                a2 = vtk.vtkActor()
                a2.SetMapper(m2)
                f.renderer.AddActor(a2)
            if first_time and self.args['axes']:
                bounds = np.array(f.clipper.bounds)
                mins = bounds[[0,2,4]]
                maxs = bounds[[1,3,5]]
                diam = np.array(maxs) - np.array(mins)
                scale = np.linalg.norm(diam, ord=-np.inf)/2
                self.oms.append(draw_axes(f.renderer, self.interactor, f.dims, scale=scale))
            if first_time:
                self.window.AddRenderer(f.renderer)

        if first_time:
            self.interactor.AddObserver('KeyPressEvent', self.key_pressed_callback)
            self.window.SetSize(2000, 1500)
            if self.args['axes']:
                for om in self.oms:
                    om.On()
                    om.InteractiveOn()
            self.interactor.Initialize()

        self.window.Render()
        if self.args['clip']:
            for f in self.frames:
                f.clipper.Activate()
        else:
            for f in self.frames:
                f.clipper.Deactivate()
        self.window.Render()

        if first_time:
            self.interactor.Start()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Visualize Slab of 4D map')
    parser.add_argument('-i', '--input', type=str, required=True, action='append', help='Name(s) of NRRD file(s) containing orbits')
    parser.add_argument('-s', '--scale', action='store_true', help='Use glyph scaling')
    parser.add_argument('-f', '--factor', type=float, default=0.01, help='Scaling factor to use')
    parser.add_argument('-m', '--min', type=int, default=2, help='Minimum orbit length')
    parser.add_argument('-c', '--color', type=str, default='index', help='Attribute to colormap (from \'index\', \'p2\', \'length\', \'random\', \'p1\', \'q1\', \'q2\', \'norm\')')
    parser.add_argument('-a', '--axes', action='store_true', help='Display coordinate axes')
    parser.add_argument('-r', '--resolution', type=int, default=0, help='Sphere resolution in theta and phi (0 to use points instead)')
    parser.add_argument('--stereo', type=float, default=None, help='Perform stereographic projection with given light position along the 4th axis')
    parser.add_argument('-x', '--max', type=int, default=None, help='Maximum number of points per orbit')
    parser.add_argument('-t', '--triangulate', action='store_true', help='Construct triangulation of point clouds')
    parser.add_argument('--coordinates', type=str, nargs=3, default=None, help='What coordinate axes to use')
    parser.add_argument('--flatten', type=float, default=None, help='Restrict fourth dimension to a thin slab')
    parser.add_argument('--clip', action='store_true', help='Toggle clipping plane')
    parser.add_argument('--total', type=int, default=None, help='Maximum total number of points to be displayed')

    args = parser.parse_args()
    print('args.input={}'.format(args.input))
    input = find_files(args.input)

    # visualize_orbit(input=input, args=vars(args))

    blah = FourDimensionalVisualizer(input, vars(args))
    blah.run()

    # visualize4d(input=input, args=vars(args))
