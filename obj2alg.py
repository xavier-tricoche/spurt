import vtk
import numpy as np

class Object2Algorithm(VTKPythonAlgorithmBase):
        def __init__(self, obj):
            if isinstance(obj, vtk.vtkPolyData):
                self.otype = 'vtkPolyData'
            elif isinstance(obj, vtk.vtkImageData) or
                 isinstance(obj, vtk.UniformGrid) or
                 isinstance(obj, vtk.vtkStructuredPoints):
                self.atype = 'vtkImageData'
            elif isinstance(obj, vtk.vtkStructuredGrid):
                self.atype = 'vtkStructuredGrid'
            elif isinstance(obj, vtk.vtkUnstructuredGrid):
                self.otype = 'vtkUnstructuredGrid'
            elif isinstance(obj, vtk.vtkPath):
                self.otype = 'vtkPath'
            else:
                raise ValueError('Unrecognized object type: {}'.format(type(obj)))

            VTKPythonAlgorithmBase.__init__(self,
                nInputPorts=0,
                nOutputPorts=1, outputType=self.otype)
            self.obj = obj

        def RequestInformation(self, request, inInfo, outInfo):
            info = outInfo.GetInformationObject(0)
            info.Set(vtkmodules.vtkCommonExecutionModel.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT(),
                (0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1), 6)
            return 1
        def RequestData(self, request, inInfo, outInfo):
            f = h5py.File("foo.h5", 'r')
            data = f['RTData'][:]
            output = dsa.WrapDataObject(vtkmodules.vtkCommonDataModel.vtkImageData.GetData(outInfo))
            output.SetDimensions(data.shape)
            output.PointData.append(data.flatten(), 'RTData')
            output.PointData.SetActiveScalars('RTData')
            return 1
