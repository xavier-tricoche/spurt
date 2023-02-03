import vtk
import numpy as np

class ClippingPlane:
    def __init__(self, plane=None):
        if plane is not None:
            self.plane = plane
        else:
            self.plane = vtk.vtkPlane()
        self.clipper = vtk.vtkClipPolyData()
        self.clipper.SetClipFunction(self.plane)
        self.clipper.InsideOutOn()
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInputConnection(self.clipper.GetOutputPort())
        self.actor = vtk.vtkActor()
        self.actor.SetMapper(self.mapper)
        self.plane_info = { 'origin': self.plane.GetOrigin(), 'normal': self.plane.GetNormal() }
        self.rep = vtk.vtkImplicitPlaneRepresentation()
        self.rep.SetPlaceFactor(1.25) # This must be set prior to placing the widget.

    def Initialize(self, input, ren, iren, **kwargs):
        print('initializing clipper in renderer {}'.format(id(ren)))
        if isinstance(input, vtk.vtkAlgorithm):
            self.clipper.SetInputConnection(input.GetOutputPort())
        else:
            self.clipper.SetInputData(input)
        ren.AddActor(self.actor)
        self.bounds = self.actor.GetBounds()
        print('actor bounds are {}'.format(self.bounds))
        self.rep.PlaceWidget(self.bounds)
        self.rep.SetNormal(self.plane.GetNormal())
        if 'origin' in kwargs.keys():
            self.rep.SetOrigin(kwargs['origin'])
        if 'normal' in kwargs.keys():
            self.rep.SetNormal(kwargs['normal'])
        if 'bounds' in kwargs.keys():
            self.bounds = kwargs['bounds']
            self.rep.PlaceWidget(self.bounds)
        if 'callback' in kwargs.keys():
            self.callback = kwargs['callback']
        else:
            self.callback = None
        # self.rep.SetRenderer(ren)

        self.widget = vtk.vtkImplicitPlaneWidget2()
        self.widget.SetInteractor(iren)
        self.widget.SetRepresentation(self.rep)
        self.widget.AddObserver('InteractionEvent', self.execute)
        self.widget.SetDefaultRenderer(ren)
        self.plane_info['origin'] = self.rep.GetOrigin()
        self.plane_info['normal'] = self.rep.GetNormal()

    def SetPositionNormal(self, position=None, normal=None):
        if position is not None:
            self.rep.SetOrigin(position)
            self.plane_info['origin'] = position
        if normal is not None:
            self.rep.SetNormal(normal)
            self.plane_info['normal'] = normal

    def Update(self, input):
        if isinstance(input, vtk.vtkAlgorithm):
            self.clipper.SetInputConnection(input.GetOutputPort())
        else:
            self.clipper.SetInputData(input)
        self.bounds = self.actor.GetBounds()

    def Deactivate(self):
        print('deactivating')
        self.plane_info['origin'] = self.plane.GetOrigin()
        self.plane_info['normal'] = self.plane.GetNormal()
        self.rep.PushPlane(+1000)
        # self.execute(self.widget, 0)
        self.plane.SetOrigin([0,0,1000])
        self.plane.SetNormal([0,0,1])

        print('plane is now set to origin={}, normal={}'.format(self.plane.GetOrigin(), self.plane.GetNormal()))
        self.Hide()
        self.widget.Off()

    def Activate(self):
        print('activating')
        self.plane.SetOrigin(self.plane_info['origin'])
        self.plane.SetNormal(self.plane_info['normal'])
        print('plane set to origin={}, normal={}'.format(self.plane.GetOrigin(), self.plane.GetNormal()))
        self.Show()
        self.widget.On()

    def Hide(self):
        self.rep.DrawPlaneOff()
        self.rep.DrawOutlineOff()

    def Show(self):
        self.rep.DrawPlaneOn()
        self.rep.DrawOutlineOn()

    def On(self):
        self.widget.On()

    def execute(self, caller, id):
        caller.GetRepresentation().GetPlane(self.plane)
        if self.callback is not None:
            self.callback(self)

if __name__ == '__main__':
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(10)
    sphere.SetThetaResolution(20)
    sphere.SetPhiResolution(20)


    clip = ClippingPlane()

    ren = vtk.vtkRenderer()
    win = vtk.vtkRenderWindow()
    win.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(win)
    clip.Initialize(sphere, ren, iren)

    win.SetSize(2000, 1000)
    iren.Initialize()
    win.Render()
    clip.On()
    iren.Start()
