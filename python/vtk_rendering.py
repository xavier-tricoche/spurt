import vtk
import numpy as np
from vtk.util import numpy_support
from numpy.random import default_rng
import os

'''
Helper function to create the mapper matching the input type.
'''

__all__ = [
    'make_mapper',
    'make_actor',
    'make_render_kit',
    'make_tubes',
    'make_spheres',
    'make_ellipsoids',
    'make_fiber_actor',
    'take_screenshot',
]

def make_mapper(input):
    if not isinstance(input, vtk.vtkAlgorithm) and not isinstance(input, vtk.vtkDataSet):
        raise ValueError(f'Invalid type ({type(intpu)}) in input of create_mapper')
    if isinstance(input, vtk.vtkAlgorithm):
        if isinstance(input, vtk.vtkPolyDataAlgorithm):
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(input.GetOutputPort())
            return mapper
        elif isinstance(input, vtk.DataSetAlgorithm):
            mapper = vtk.vtkDataSetMapper()
            mapper.SetInputConnection(input.GetOutputPort())
            return mapper
        else:
            raise ValueError(f'Unrecognized algorithm type ({type(input)}) in create_mapper')
    elif isinstance(input, vtk.vtkPolyData):
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(input)
        return mapper
    elif isinstance(input, vtk.vtkDataSet):
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(input)
        return mapper
    else:
        raise ValueError(f'Unrecognized dataset type ({type(input)}) in create_mapper')

'''
Helper function to create the actor matching the input type.
'''
def make_actor(data, ctf=None, show_scalars=False, **kwargs):
    if isinstance(data, vtk.vtkMapper):
        mapper = data
    else: mapper = vtk.vtkDataSetMapper()
    if ctf is not None:
        mapper.SetLookupTable(ctf)
        mapper.ScalarVisibilityOn()
    elif show_scalars:
        mapper.ScalarVisibilityOn()
    else:
        mapper.ScalarVisibilityOff()
    if isinstance(data, vtk.vtkAlgorithm):
        mapper.SetInputConnection(data.GetOutputPort())
    else:
        mapper.SetInputData(data)
    prop = vtk.vtkProperty(**kwargs)
    return vtk.vtkActor(mapper=mapper, property=prop)

'''
Helper function to create the rendering elements (renderer, window, interactor).
'''
def make_render_kit(actors=[], **kwargs):
    renderer = vtk.vtkRenderer(background=kwargs.get('background', [0,0,0]))
    for a in actors:
        if isinstance(a, vtk.vtkActor2D):
            renderer.AddActor2D(a)
        else:
            renderer.AddActor(a)
    #window = vtk.vtkRenderWindow(size=kwargs.get('size', (1024, 1024)))
    window = vtk.vtkRenderWindow(size=(1024, 1024))
    window.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor(render_window=window)
    return renderer, window, interactor

'''
Helper function to wrap source geometry with tubes
'''
def make_tubes(source, radius=0.1, resolution=12):
    if not isinstance(source, vtk.vtkAlgorithm):
        s = vtk.vtkTrivialProducer()
        s.SetOutput(source)
        source = s
    tubes = vtk.vtkTubeFilter(radius=radius, number_of_sides=resolution)
    tubes.SetInputConnection(source.GetOutputPort())
    return make_actor(tubes)

'''
Helper function to create spherical glyphs
'''
def make_spheres(source, radius=0.2, resolution=12):
    if not isinstance(source, vtk.vtkAlgorithm):
        data = source
        source = vtk.vtkTrivialProducer()
        source.SetOutput(data)
    else:
        data = source.GetOutputDataObject()
    if data.GetVerts() is None:
        cells = vtk.vtkCellArray()
        n = data.GetNumberOfPoints()
        for i in range(n):
            data.InsertNextCell(1)
            data.InsertCellPoint(i)
        data.SetVerts(cells)
    sphere = vtk.vtkSphereSource(theta_resolution=resolution, phi_resolution=resolution, radius=radius)
    glyph = vtk.vtkGlyph3D(scaling=False)
    glyph.SetSourceConnection(sphere.GetOutputPort())
    glyph.SetInputConnection(source.GetOutputPort())
    return make_actor(glyph, True)

'''
Make tensor ellipsoid glyphs
'''
def make_ellipsoids(source, scaling=10, resolution=18):
    if not isinstance(source, vtk.vtkAlgorithm):
        data = source
        source = vtk.vtkTrivialProducer()
        source.SetOutput(data)
    else:
        data = source.GetOutputDataObject()
    sphere = vtk.vtkSphereSource(theta_resolution=resolution, phi_resolution=resolution)
    tglyph = vtk.vtkTensorGlyph(scale_factor=scaling, clamp_scaling=True, three_glyphs=False)
    tglyph.SetSourceConnection(sphere.GetOutputPort())
    tglyph.SetInputConnection(source.GetOutputPort())
    return make_actor(tglyph, True)

'''
Make actor from single fiber
'''
def make_fiber_actor(fiber, values=None, radius=0.1, resolution=12, as_tube=True):
    poly = make_points(fiber)
    if values is not None:
        if not isinstance(values, vtk.vtkObject):
            values = numpy_support.numpy_to_vtk(values)
        poly.GetPointData().SetScalars(values)
    cells = vtk.vtkCellArray()
    cells.InsertNextCell(len(fiber))
    for i in range(len(fiber)):
        cells.InsertCellPoint(i)
    poly.SetLines(cells)
    if as_tube:
        a, tubes = make_tubes(poly, radius=radius, resolution=resolution)
    else:
        m = vtk.vtkPolyDataMapper(scalar_visibility=True)
        m.SetInputData(poly)
        a = vtk.vtkActor(mapper=m)
        a.GetProperty().RenderLinesAsTubesOn()
        a.GetProperty().SetLineWidth(radius)
        tubes = poly
    if values is None:
        rng = default_rng()
        col = rng.random([3])
        a.GetProperty().SetColor(col)
    else:
        a.GetMapper().ScalarVisibilityOn()
    return a, tubes

def take_screenshot(window, filename, format=None):
    if format is None and os.path.splitext(filename)[1]:
        format = os.path.splitext(filename)[1][1:]
    elif format is None:
        format = 'jpg'
        filename += '.jpg'
    image = vtk.vtkWindowToImageFilter()
    image.SetInput(window)
    match format:
        case 'png':
            writer = vtk.vtkPNGWriter(file_name=filename)
        case 'jpg' | 'jpeg':
            writer = vtk.vtkJPEGWriter(file_name=filename)
        case 'tif' | 'tiff':
            writer = vtk.vtkTIFFWriter(file_name=filename)
        case 'bmp':
            writer = vtk.vtkBMPWriter(file_name=filename)
        case _:
            raise ValueError(f"Unsupported format: {format}")
    writer.SetInputConnection(image.GetOutputPort())
    window.Render()
    writer.Write()

def main():
    src = vtk.vtkSphereSource(theta_resolution=100, phi_resolution=100)
    actor = make_actor(src)
    actor.GetProperty().SetColor(1,0,0)

    src.Update()
    copy = vtk.vtkPolyData()
    copy.DeepCopy(src.GetOutput())
    xform = vtk.vtkTransform()
    xform.Translate(1,1,0)
    xformpd = vtk.vtkTransformPolyDataFilter()
    xformpd.SetInputData(copy)
    xformpd.SetTransform(xform)
    xformpd.Update()
    actor2 = make_actor(xformpd.GetOutput())
    actor2.GetProperty().SetColor(0,0,1)

    image = vtk.vtkImageData(dimensions=[1000, 1000, 1], origin=[0,0,0], spacing=[0.001, 0.001, 1])
    t = np.linspace(0, 2 * np.pi, 1000)
    data2d = np.sin(t)[:, np.newaxis] * np.cos(t)[np.newaxis, :]
    image.GetPointData().SetScalars(numpy_support.numpy_to_vtk(data2d.flatten('F')))
    actor3 = make_actor(image)

    renderer = vtk.vtkRenderer(background=[0,0,0])
    renderer.AddActor(actor)
    renderer.AddActor(actor2)
    renderer.AddActor(actor3)

    window = vtk.vtkRenderWindow(size=[1920,1080])
    window.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor(render_window=window)
    interactor.Initialize()
    window.Render()
    interactor.Start()

if __name__ == '__main__':
    main()