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
import vtk_colors

from vtk.util.numpy_support import *

with_axes = False

name_to_column = {
    'p1': 'x',
    'p2': 'y',
    'q1': 'xdot',
    'q2': 'ydot',
    'r' : 'norm'
}


def draw_axes(renderer, dims):
    arrows = []
    for i in range(3):
        a = vtk.vtkArrowSource()
        a.SetTipResolution(20)
        a.SetShaftResolution(10)
        a.SetTipRadius(0.01)
        a.SetTipLength(0.035)
        a.SetShaftRadius(0.002)
        arrows.append(a)

    mx = vtk.vtkPolyDataMapper()
    ax = vtk.vtkActor()
    ax.SetMapper(mx)
    mx.SetInputConnection(arrows[0].GetOutputPort())
    ax.GetProperty().SetColor(1,0,0)

    rz = vtk.vtkTransform()
    rz.Identity()
    rz.RotateZ(90)
    ry = vtk.vtkTransform()
    ry.Identity()
    ry.RotateY(-90)

    ty = vtk.vtkTransformFilter()
    ty.SetTransform(rz)
    ty.SetInputConnection(arrows[1].GetOutputPort())
    my = vtk.vtkPolyDataMapper()
    ay = vtk.vtkActor()
    ay.SetMapper(my)
    my.SetInputConnection(ty.GetOutputPort())
    ay.GetProperty().SetColor(0,1,0)

    tz = vtk.vtkTransformFilter()
    tz.SetTransform(ry)
    tz.SetInputConnection(arrows[2].GetOutputPort())
    mz = vtk.vtkPolyDataMapper()
    az = vtk.vtkActor()
    az.SetMapper(mz)
    mz.SetInputConnection(tz.GetOutputPort())
    az.GetProperty().SetColor(0,0,1)

    renderer.AddActor(ax)
    renderer.AddActor(ay)
    renderer.AddActor(az)

def mindist(points, nneighs=2):
    kdtree = spatial.KDTree(points)
    dist, ids = kdtree.query(points, k=nneighs+1)
    # print('dist={}'.format(dist))
    # print('ids={}'.format(ids))
    dist = dist[:,1:]
    # print('dist={}'.format(dist))
    max = np.max(dist)
    print('max 2nd min distance is {}'.format(max))
    return max

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
        # glyphs.SetScaleModeToScaleByScalar()
        glyphs.SetColorModeToColorByScalar()
        # glyphs.SetScaleFactor(scale_factor)
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
        for fname in filenames:
            fname = fname.name
            print(f'fname={fname}')
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

def visualize_orbit(input, scale=False, scale_factor=0.001, colorby=None, sres=6):
    renderer = vtk.vtkRenderer()
    window = vtk.vtkRenderWindow()
    window.AddRenderer(renderer)
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
    for i, name in enumerate(names):
        # print('processing {}'.format(name))
        data, header = nrrd.read(name)
        df = pd.DataFrame(data, index=['x', 'y', 'xdot', 'ydot']).transpose()
        df['order'] = np.linspace(0, 1, data.shape[1])
        if colorby == 'random':
            df['index'] = i
        if main_df is not None:
            main_df = pd.concat([main_df, df], ignore_index=True)
        else:
            main_df = df
    main_df.reset_index()

    if colorby == 'r' or colorby == 'norm':
         main_df['norm'] = main_df.apply(np.linalg.norm, axis=1)

    c = numpy_to_vtk(main_df[['x', 'y', 'xdot']].to_numpy())
    points = vtk.vtkPoints()
    points.SetData(c)
    cloud = vtk.vtkPolyData()
    cloud.SetPoints(points)
    if colorby in ['p1', 'p2', 'q1', 'q2', 'r']:
        attribute = numpy_to_vtk(main_df[name_to_column[colorby]].to_numpy())
        attribute.SetName(name_to_column[colorby])
    elif colorby in name_to_column.keys():
        attribute = numpy_to_vtk(main_df[colorby].to_numpy())
        attribute.SetName(colorby)
    elif colorby in ['r', 'norm']:
        attribute = numpy_to_vtk(main_df['norm'].to_numpy())
        attribute.SetName('norm')
    elif colorby in ['index', 'order']:
        attribute = numpy_to_vtk(main_df['order'].to_numpy())
        attribute.SetName('order')
    elif colorby == 'random':
        attribute = numpy_to_vtk(main_df['index'])
    cloud.GetPointData().AddArray(attribute)
    cloud.GetPointData().SetActiveScalars(attribute.GetName())
    glyphs = make_glyphs(cloud, scale, scale_factor, sres)

    m = vtk.vtkPolyDataMapper()
    a = vtk.vtkActor()
    a.SetMapper(m)
    a.GetProperty().RenderPointsAsSpheresOn()
    a.GetProperty().SetPointSize(10)

    m.SetInputConnection(glyphs.GetOutputPort())

    print(f'the range of values in selected attribute is {attribute.GetRange()}')

    if colorby == 'order' or colorby == 'index' or colorby == 'random':
        m.ScalarVisibilityOn()
        cmap = vtk_colors.make_colormap('viridis', attribute.GetRange())
        m.SetLookupTable(cmap)
    elif colorby in ['p1', 'p2', 'q1', 'q2', 'r'] or colorby in name_to_column.keys():
        m.ScalarVisibilityOn()
        range = attribute.GetRange()
        if colorby not in ['r', 'norm']:
            cmap = vtk.vtkColorTransferFunction()
            cmap.AddRGBPoint(range[0], 0, 0, 1)
            cmap.AddRGBPoint(0, 1, 1, 1)
            cmap.AddRGBPoint(range[1], 1, 0, 0)
        else:
            cmap = vtk_colors.make_colormap('viridis', [range[0], range[1]])
        m.SetLookupTable(cmap)
    else:
        m.ScalarVisibilityOff()
        a.GetProperty().SetColor(1,1,0)
    renderer.AddActor(a)

    global with_axes
    if with_axes:
        draw_axes(renderer, ['x', 'y', 'xdot'])

    window.SetSize(2000, 1500)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)
    interactor.Initialize()
    window.Render()
    interactor.Start()


def visualize_4d(name, scale=False, scale_factor=0.001, min_length=2, colorby=None, sres=6):
    # x, y, xdot, ydot
    # x, y, xd
    # x, y, yd
    # x, xd, yd
    # y, xd, yd
    data, header = nrrd.read(name)
    df = pd.DataFrame(data, index = ['x', 'y', 'xdot', 'ydot', 'index']).transpose()
    indices = df['index'].unique()

    random.seed(a=13081975)
    unique_colors = np.array([ random.randrange(0,255) for i in range(3*len(indices)) ])
    print(f'There are {len(indices)} unique indices')
    colors = numpy.array([ [ unique_colors[int(3*i)], unique_colors[int(3*i+1)], unique_colors[int(3*i+2)] ] for i in df['index'] ])
    print('colors={}'.format(colors))
    # df['colors'] = colors

    array = colors
    z = np.asarray(array)
    print('z={}'.format(z))
    if not z.flags.contiguous:
        z = np.ascontiguousarray(z)

    shape = z.shape
    assert z.flags.contiguous, 'Only contiguous arrays are supported.'
    assert len(shape) < 3, \
           "Only arrays of dimensionality 2 or lower are allowed!"
    assert not numpy.issubdtype(z.dtype, numpy.dtype(complex).type), \
           "Complex numpy arrays cannot be converted to vtk arrays."\
           "Use real() or imag() to get a component of the array before"\
           " passing it to vtk."

    print('the type is {}'.format(z.dtype))

    # # First create an array of the right type by using the typecode.
    # if array_type:
    #     vtk_typecode = array_type
    # else:
    #     vtk_typecode = get_vtk_array_type(z.dtype)
    # result_array = create_vtk_array(vtk_typecode)

    print(array)

    npts = data.shape[1]
    ids = numpy_to_vtk(df['index'].to_numpy())
    dist = numpy_to_vtk(df['ydot'].to_numpy())
    colors = numpy_to_vtk(colors, array_type=vtk.VTK_UNSIGNED_CHAR)

    coherence = []
    for i in indices:
        subdf = df[df['index']==i]
        coherence.append(mindist(subdf[['x', 'y', 'xdot', 'ydot']]))
        print('coherence = {}'.format(coherence))

    # ids = vtk.vtkFloatArray()
    # ids.SetNumberOfTuples(npts)
    # ids.SetNumberOfComponents(1)
    # ids.SetName('id')
    #
    # dist = vtk.vtkFloatArray()
    # dist.SetNumberOfTuples(npts)
    # dist.SetNumberOfComponents(1)
    # dist.SetName('p2')
    #
    # counter = 0
    # colors = vtk.vtkUnsignedCharArray()
    # vtk_colors.SetNumberOfComponents(3)
    # vtk_colors.SetNumberOfTuples(npts)
    # vtk_colors.SetName('colors')


    # coherence = []
    #
    # for k in valid.keys():
    #     sub = df[df['index']==k]
    #     col = [random.randrange(0,255), random.randrange(0,255), random.randrange(0,255)]
    #     for i, row in tqdm(sub.iterrows()):
    #         ids.SetTuple1(counter, row['index'])
    #         dist.SetTuple1(counter, row['ydot'])
    #         vtk_colors.SetTuple3(counter, col[0], col[1], col[2])
    #         counter = counter+1
    #     coherence.append(mindist(sub[['x', 'y', 'xdot', 'ydot']))

    coherence_ranking = np.argsort(coherence)

    clouds = []
    all_dims = [ ['x', 'y', 'xdot'], ['x', 'y', 'ydot'],
             ['x', 'xdot', 'ydot'], ['y', 'xdot', 'ydot'] ]

    xmins = [ 0, 0.5, 0, 0.5 ]
    ymins = [ 0, 0, 0.5, 0.5 ]

    window = vtk.vtkRenderWindow()

    for n in range(4):
        dims = all_dims[n]
        c = numpy_to_vtk(df[[dims[0], dims[1], dims[2]]].to_numpy())
        # c = vtk.vtkFloatArray()
        # c.SetNumberOfComponents(3)
        # c.SetNumberOfTuples(npts)
        # counter = 0
        # dims = all_dims[n]
        # for k in valid.keys():
        #     sub = df[df['index']==k]
        #     for i, row in tqdm(sub.iterrows()):
        #         c.SetTuple3(counter, row[dims[0]], row[dims[1]], row[dims[2]])
        #         counter = counter+1

        points = vtk.vtkPoints()
        points.SetData(c)
        cloud = vtk.vtkPolyData()
        cloud.SetPoints(points)
        cloud.GetPointData().AddArray(ids)
        cloud.GetPointData().AddArray(dist)
        cloud.GetPointData().AddArray(colors)
        cloud.GetPointData().SetActiveScalars('colors')
        clouds.append(cloud)
        glyphs = make_glyphs(cloud, scale, scale_factor, sres)

        renderer = vtk.vtkRenderer()
        window.AddRenderer(renderer)
        print(f'n={n}')
        renderer.SetViewport(xmins[n], ymins[n], xmins[n]+0.5, ymins[n]+0.5)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(glyphs.GetOutputPort())
        mapper.ScalarVisibilityOn()

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().VertexVisibilityOn()
        renderer.AddActor(actor)
        draw_axes(renderer, dims)


    window.SetSize(2000, 1500)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)
    interactor.Initialize()
    window.Render()
    interactor.Start()


def visualize(input, scale=False, scale_factor=0.001,
              min_length=2, colorby=None, sres=6):


    hits_name = None
    orbits_name = None
    if input.find('hits') >= 0:
        orbits_name = input.replace('hits', 'orbits')
        hits_name = input
    elif input.find('orbits') >= 0:
        hits_name = input.replace('orbits', 'hits')
        orbits_name = input
    elif input.find('XX'):
        hits_name = input.replace('XX', 'hits')
        orbits_name = input.replace('XX', 'orbits')
    else:
        raise ValueError('Invalid Input Filename: {}'.format(input))

    visualize_4d(orbits_name, scale, scale_factor, min_length, colorby, sres)
    return 0



    hits_data, header = nrrd.read(hits_name)
    print(hits_data.shape)
    # print(data)
    # print(header)
    nhits= hits_data.shape[1]

    # figure out which orbits are valid
    hits_df = pd.DataFrame(hits_data, index = ['x', 'y', 'xdot', 'ydot', 'index']).transpose()
    print(df)
    indices = hits_df['index'].unique()
    print(f'There are {len(indices)} unique indices')
    lengths = [ 0 ] * len(indices)
    valid = {}
    actual_nhits = 0
    for i in indices:
        subdf = hits_df[df['index']==i]
        n = subdf.shape[0]
        lengths[int(i)] = n
        if n >= min_length:
            valid[i] = n
            actual_nhits += n

    coords = vtk.vtkFloatArray()
    coords.SetNumberOfComponents(3)
    coords.SetNumberOfTuples(actual_nhits)

    ids = vtk.vtkFloatArray()
    ids.SetNumberOfTuples(actual_nhits)
    ids.SetNumberOfComponents(1)
    ids.SetName('id')

    dist = vtk.vtkFloatArray()
    dist.SetNumberOfTuples(actual_nhits)
    dist.SetNumberOfComponents(1)
    dist.SetName('p2')

    leng = vtk.vtkFloatArray()
    leng.SetNumberOfTuples(actual_nhits)
    leng.SetNumberOfComponents(1)
    leng.SetName('length')

    counter = 0
    for k in valid.keys():
        sub = hits_df[hits_df['index']==k]
        for i, row in sub.iterrows():
            coords.SetTuple3(counter, row['x'], row['y'], row['xdot'])
            # print(f'set #{counter} to {row['x']}, {row['y']}, {row['xdot']}')
            ids.SetTuple1(counter, row['index'])
            dist.SetTuple1(counter, row['ydot'])
            leng.SetTuple1(counter, valid[k])
            counter = counter+1
        mindist(sub[['x', 'y', 'xdot']])

    points = vtk.vtkPoints()
    points.SetData(coords)
    geom = vtk.vtkPolyData()
    geom.SetPoints(points)
    geom.GetPointData().AddArray(ids)
    geom.GetPointData().AddArray(dist)
    geom.GetPointData().AddArray(leng)
    if colorby is None:
        colorby = 'length'
    geom.GetPointData().SetActiveScalars(colorby)
    color_range = geom.GetPointData().GetScalars().GetRange()
    print(f'color_range={color_range}')

    geom.ComputeBounds()

    print(f'the bounds of the polydata are {geom.GetBounds()}')

    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(0.001)
    sphere.SetThetaResolution(sres)
    sphere.SetPhiResolution(sres)
    glyphs = vtk.vtkGlyph3D()
    glyphs.SetSourceConnection(sphere.GetOutputPort())
    glyphs.SetInputData(geom)
    glyphs.ScalingOn()
    glyphs.SetScaleModeToScaleByScalar()
    glyphs.SetColorModeToColorByScalar()
    glyphs.SetScaleFactor(scale_factor)

    ctf = vtk.vtkColorTransferFunction()
    ctf.AddRGBPoint(color_range[0], 0, 0, 1)
    ctf.AddRGBPoint(color_range[1], 1, 1, 0)
    # print(f'color ranges from {np.min(df[colorby])} and {np.max(df[colorby])}')
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyphs.GetOutputPort())
    mapper.ScalarVisibilityOn()
    mapper.SetLookupTable(ctf)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    cube = vtk.vtkCubeSource()
    cube.SetXLength(1)
    cube.SetYLength(1)
    cube.SetZLength(1)
    cube.SetCenter(0,0,0)

    cube_mapper = vtk.vtkPolyDataMapper()
    cube_mapper.SetInputConnection(cube.GetOutputPort())
    cube_actor = vtk.vtkActor()
    cube_actor.SetMapper(cube_mapper)
    cube_actor.GetProperty().SetColor(1,1,0)
    cube_actor.GetProperty().SetRepresentationToWireframe()
    cube_actor.GetProperty().RenderLinesAsTubesOn()
    # cube_actor.GetProperty().SetLineWidth(0.5)

    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.AddActor(cube_actor)
    renderer.ResetCamera()

    window = vtk.vtkRenderWindow()
    window.AddRenderer(renderer)
    window.SetSize(1280,900)

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)
    interactor.Initialize()

    window.Render()
    interactor.Start()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Visualize Slab of 4D map')
    parser.add_argument('-i', '--input', type=str, required=True, action='append', help='Name(s) of NRRD file(s) containing orbits')
    parser.add_argument('-s', '--scale', action='store_true', help='Use glyph scaling')
    parser.add_argument('-f', '--factor', type=float, default=0.01, help='Scaling factor to use')
    parser.add_argument('-m', '--min', type=int, default=2, help='Minimum orbit length')
    parser.add_argument('-c', '--color', type=str, default='index', help='Attribute to colormap (from \'index\', \'p2\', \'length\', \'random\', \'p1\', \'q1\', \'q2\', \'norm\')')
    parser.add_argument('-a', '--axes', action='store_true', help='Display coordinate axes')
    parser.add_argument('-r', '--resolution', type=int, default=6, help='Sphere resolution in theta and phi')

    args = parser.parse_args()
    print('args.input={}'.format(args.input))
    input = find_files(args.input)
    # print('we found\n\n {} \nfrom\n\n {}'.format(input, args.input))

    # print(f'input={input}')
    with_axes = args.axes

    # visualize_orbit(input=args.input, scale=args.scale, scale_factor=args.factor, min_length=args.min, colorby = args.color, sres=args.resolution)
    visualize_orbit(input=input, scale=args.scale, scale_factor=args.factor, colorby = args.color, sres=args.resolution)
