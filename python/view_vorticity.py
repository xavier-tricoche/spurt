import numpy as np 
import vtk 
from vtk.util import numpy_support as nps
import nrrd 
import vtk_colors_updated as colors
import vtk_colorbar_updated as bar
import vtk_rendering as ren

w = 7.292123516990375e-05 # Earth sideral spin rate 

def main():
    import argparse 
    parser = argparse.ArgumentParser(description='Color map vorticity')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input filename')
    parser.add_argument('-o', '--output', type=str, help='Output filename')
    parser.add_argument('-r', '--range', type=float, nargs=2, help='Bounds of color map')
    parser.add_argument('-c', '--coriolis', action='store_true', help='Apply Coriolis term')
    parser.add_argument('-l', '--lat', type=float, help='Latitude of reference point')
    args = parser.parse_args()

    data, header = nrrd.read(args.input)

    if args.coriolis and args.lat is not None:
        lat = np.deg2rad(args.lat)
        data = data / (2*w*np.sin(lat))

    if args.range is not None:
        minval = args.range[0]
        maxval = args.range[1] 
    else:
        minval = np.min(data)
        maxval = np.max(data)

    img = vtk.vtkImageData()
    img.SetDimensions(data.shape[0], data.shape[1], 1)
    #img.SetSpacing(header['spacings'][0], header['spacings'][1], 1.)
    img.SetSpacing(0.7, 1, 1)
    #img.SetOrigin(header['axis mins'][0], header['axis mins'][1], 0.)
    img.GetPointData().SetScalars(nps.numpy_to_vtk(data.ravel('F')))

    ctf = colors.make_colormap('seismic', [minval, 0, maxval], diverging=True)
    actor = ren.make_actor(img, ctf)
    params = bar.ColorbarParam(title='w/f', title_col=[0,0,0], title_font_size=30, 
        label_col=[0,0,0], pos=[0.95, 0.75], width=50, height=300, 
        nlabels=5, font_size=22, title_offset=10, )
    colorbar = bar.Colorbar(ctf, param=params)

    renderer, window, interactor = ren.make_render_kit(actors=[actor, colorbar.get()], size=data.shape)
    renderer.SetBackground(1, 1, 1)
    interactor.Initialize()
    renderer.ResetCamera()
    window.Render()
    interactor.Start()
    
        
if __name__ == '__main__':
    main()