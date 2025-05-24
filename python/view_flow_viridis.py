import numpy as np
import vtk
from vtk.util import numpy_support as nps
import argparse 
import sys
import nrrd
import time
import os
import vtk_camera
import vtk_colors
import regex as re

from scipy.interpolate import RegularGridInterpolator

def make2dcolors(width=64, height=64, color1=[0,0,255], color2=[255,0,0]):
    init = np.zeros((2,2,3))
    init[0,0,:] = [0,0,0]
    init[1,0,:] = color1 
    init[0,1,:] = color2
    init[1,1,:] = [255, 255, 255]

    interp = RegularGridInterpolator(([0., 1.], [0., 1.]), init)

    x = np.linspace(start=0, stop=1, num=width, endpoint=True)
    y = np.linspace(start=0, stop=1, num=height, endpoint=True)
    pts = [np.array([u,v]) for v in y for u in x]
    r = interp(pts)
    data = np.array(r).astype(np.uint8)
    return data

def key_pressed_callback(obj, event):
    # ---------------------------------------------------------------
    # Attach actions to specific keys
    # ---------------------------------------------------------------
    key = obj.GetKeySym()
    if key == "h":
        print("Commands:\n 's': save frame\n 'c': save camera setting\n 'h': print this message\n 'q': quit the program")
    elif key == "s":
        save_frame()
    elif key == "c":
        vtk_camera.save_camera(renderer=obj.GetRenderWindow().GetRenderers().GetFirstRenderer())
        print(f'Camera settings saved to camera.json')
    elif key == "q":
        if args.verbose:
            print("User requested exit.")
        sys.exit()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Visualize path of particles')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input filename')
    parser.add_argument('-r', '--res', type=int, nargs=2, default=[2000, 1600], help='Window resolution')
    parser.add_argument('-p', '--path', type=str, help='Path of data files to prepend to their name')
    parser.add_argument('-o', '--output', type=str, help='Output basename for rendered frames')
    parser.add_argument('-c', '--color', type=str, default='value', help='Attribute to use for color coding')
    parser.add_argument('-n', '--normalize', action='store_true', help='Apply time normalization to LAVD values')
    parser.add_argument('--camera', type=str, help='Input camera settings')
    parser.add_argument('-v', '--verbose', action='store_true', help='Toggle verbose mode')

    args = parser.parse_args()

    if os.path.splitext(args.input)[1] == '.txt':
        print(f'import filenames from {args.input}')
        step_names = []
        with open(args.input) as fp:
            step_names = fp.readlines()
    else:
        step_names = [ args.input ]

    pos = []
    vals = []
    names = []
    maxs = []
    times = []

    make2dcolors()

    verts = vtk.vtkCellArray()
    
    single_frame = len(step_names) == 1
    timestr = re.compile(r".*([0]*[0-9]{5})h")

    for i, name in enumerate(step_names):
        if name[-1] == '\n':
            name = name[:-1]
        times.append(int(timestr.search(name).group(1)))
        
        if not os.path.exists(name) and os.path.split(name)[0] == '':
            if args.path is not None:
                name = os.path.join(args.path, name)
            elif os.path.split(args.input)[0] != '':
                name = os.path.join(os.path.split(args.input)[0], name)
        if args.verbose: print(f'name={name}')
        names.append(name)
        data, _ = nrrd.read(name) # [x, y, lavd] x M x N
        vals.append(np.copy(data[2,:,:].ravel(order='F')))
        vals[-1] /= float(times[-1])
        data[2,:,:] = 0
        if not i:
            npts = data.shape[1]*data.shape[2]
            width = data.shape[1]
            height = data.shape[2]
            for n in range(npts):
                verts.InsertNextCell(1)
                verts.InsertCellPoint(n)
        pos.append(np.swapaxes(data.reshape((3,-1), order='F'), 0, 1))
        maxs.append(np.max(vals[-1]))
        if args.verbose:
            print(f'Step #{i}, name: {name}, bounding box: [{np.min(pos[-1][0,:])}, {np.min(pos[-1][1,:])}] -> [{np.max(pos[-1][0,:])}, {np.max(pos[-1][1,:])}], value range: [{np.min(vals[-1])}, {np.max(vals[-1])}')

    allvalues = np.concatenate(vals, axis=0)
    print(allvalues.shape)
    invalid = np.less(allvalues, 0)
    valid = np.delete(allvalues, invalid)
    print(f'there are {np.nonzero(invalid)[0].shape} invalid values')
    print(f'after removing them, there are {valid.shape} valid values left')
    print(f'percentile of all valid values: {np.percentile(valid, [0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100])}')
    ctrl_pts = np.percentile(valid, [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])


    maxvalue = np.max(maxs)
    if args.color == 'value':
        ctf = vtk_colors.make_colormap('viridis', ctrl_pts)
        ctf.AddRGBPoint(0, 0, 0, 0)
    elif args.color == 'pos' or args.color == 'position':
        colors = make2dcolors(width, height, color1=[0,0,255], color2=[255, 255, 0])

    if args.color != 'value':
        print(colors.shape)
        colors = nps.numpy_to_vtk(colors)
        
    ren = vtk.vtkRenderer()
    win = vtk.vtkRenderWindow()
    win.AddRenderer(ren)
    ren.SetBackground(0, 0, 0)
    win.SetSize(args.res[0], args.res[1])
    mapper = vtk.vtkPolyDataMapper()
    mapper.ScalarVisibilityOn()
    if args.color == 'value':
        mapper.SetLookupTable(ctf)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetPointSize(2)
    actor.GetProperty().RenderPointsAsSpheresOn()
    ren.AddActor(actor)
    if args.camera is not None:
        ren.SetActiveCamera(vtk_camera.load_camera(args.camera))

    if args.output is not None:
        capture = vtk.vtkWindowToImageFilter()
        capture.SetInput(win)
        writer = vtk.vtkJPEGWriter()
        writer.SetInputConnection(capture.GetOutputPort())

    for p, f, name in zip(pos, vals, names):
        # print(f'p={p},\n {p.shape}')
        if args.verbose: print(name)
        shortname = os.path.splitext(os.path.split(name)[1])[0]
        coords = nps.numpy_to_vtk(p)
        pts = vtk.vtkPoints()
        pts.SetData(coords)
        scl = nps.numpy_to_vtk(f)
        poly = vtk.vtkPolyData()
        poly.SetPoints(pts)
        poly.GetPointData().SetScalars(scl)
        poly.SetVerts(verts)
        if args.color != 'value':
            poly.GetPointData().SetScalars(colors)
        mapper.SetInputData(poly)
        if args.camera is None:
            ren.ResetCamera()
        if single_frame:
            inter = vtk.vtkRenderWindowInteractor()
            inter.SetRenderWindow(win)
            inter.Initialize()
        win.Render()
        if single_frame and args.output is None:
            inter.AddObserver("KeyPressEvent", key_pressed_callback)
            inter.Start()
        elif args.output is not None:
            win.Render()
            capture.Update()
            writer.SetFileName(f'{args.output}_{shortname}.jpg')
            writer.Write()
        else:
            time.sleep(1)



    


    
        


