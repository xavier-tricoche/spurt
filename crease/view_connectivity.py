import vtk 
import numpy as np 
from vtk.util import numpy_support as nps 
import argparse 
import sys 
import os 
import networkx as nx 
from tqdm import tqdm
import pandas as pd
import pickle

def cell_neighbors(df, cellid):
    pts = df[df['cellid'] == cellid]
    ids = [0, 0, 0]
    ids[0] = pts.iloc[0, 1]
    ids[1] = pts.iloc[0, 2]
    ids[2] = pts.iloc[0, 3]

    selection = df[ (df['pt0'] == ids[0]) | 
                    (df['pt1'] == ids[0]) | 
                    (df['pt2'] == ids[0]) |
                    (df['pt0'] == ids[1]) | 
                    (df['pt1'] == ids[1]) | 
                    (df['pt2'] == ids[1]) |
                    (df['pt0'] == ids[2]) | 
                    (df['pt1'] == ids[2]) | 
                    (df['pt2'] == ids[2])
                  ]
    selection = selection[selection['cellid'] != cellid]
    array = selection.to_numpy()
    match = (array[:,1:] == ids[0]) | (array[:,1:] == ids[1]) | (array[:,1:] == ids[2]).astype(int)
    match = np.sum(match, axis=1)
    return [ selection.iloc[i,0] for i in range(array.shape[0]) if match[i] == 2 ]

def compute_connectivity_graph(data):
    cellids = np.zeros(data.GetNumberOfCells(), dtype=int)
    pt0 = np.zeros(data.GetNumberOfCells(), dtype=int)
    pt1 = np.zeros(data.GetNumberOfCells(), dtype=int)
    pt2 = np.zeros(data.GetNumberOfCells(), dtype=int)

    ids = vtk.vtkIdList()
    for n in tqdm(range(cellids.shape[0])):
        data.GetCellPoints(n, ids)
        cellids[n] = n 
        pt0[n] = ids.GetId(0)
        pt1[n] = ids.GetId(1)
        pt2[n] = ids.GetId(2)

    df = pd.DataFrame({'cellid': cellids, 'pt0': pt0, 'pt1': pt1, 'pt2': pt2})

    seen = np.zeros(data.GetNumberOfCells(), dtype=bool)
    G = nx.Graph()
    for i in tqdm(range(seen.shape[0])):
        seen[i] = True
        neighs = cell_neighbors(df, i)
        for neigh in neighs:
            G.add_edge(i, neigh)
    return G

def get_connectivity_graph(data, filename):
    if os.path.isfile(filename):
        G = pickle.load(open(filename, 'rb'))
    else:
        G = compute_connectivity_graph(data)
        pickle.dump(G, open(filename, 'wb'))
    return G

def get_actor(data):
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(data)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    return actor

# Catch mouse events
class MouseInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):
    def __init__(self, data, filename, renderer):
        # self.AddObserver('LeftButtonPressEvent', self.left_button_press_event)
        self.AddObserver('RightButtonPressEvent', self.right_button_press_event)
        self.data = data
        self.src_poly = vtk.vtkPolyData()
        self.src_actor = get_actor(self.src_poly)
        self.src_actor.GetProperty().SetColor(1,0,0)
        self.src_actor.VisibilityOff()
        self.renderer = renderer
        renderer.AddActor(self.src_actor)
        self.trg_poly = vtk.vtkPolyData()
        self.trg_actor = get_actor(self.trg_poly)
        self.trg_actor.GetProperty().SetColor(0,1,0)
        self.trg_actor.VisibilityOff()
        renderer.AddActor(self.trg_actor)
        self.path_poly = vtk.vtkPolyData()
        self.path_actor = get_actor(self.path_poly)
        self.path_actor.GetProperty().SetColor(1,0,1)
        self.path_actor.GetProperty().SetLineWidth(4)
        self.path_actor.VisibilityOff()
        renderer.AddActor(self.path_actor)
        self.G = get_connectivity_graph(self.data, filename)
        self.source = -1
        self.target = -1
        self.picker = vtk.vtkCellPicker()
        self.picker.SetTolerance(0.0005)

    def anyevent(self, obj, event):
        print(f'event is {event}')
        self.OnMouseEvent()

    def lift_coords(self, p):
        p0 = np.array(self.renderer.GetActiveCamera().GetPosition())
        q = np.array(p)
        d = p0 - q
        d/=np.linalg.norm(d)
        return q + 0.001 * d

    def draw_path(self, path):
        cells = vtk.vtkCellArray()
        coords = vtk.vtkFloatArray()
        coords.SetNumberOfComponents(3)
        coords.SetNumberOfTuples(len(path))
        ids = vtk.vtkIdList()
        for i, cellid in enumerate(path):
            self.data.GetCellPoints(cellid, ids)
            p = np.zeros(3, dtype=float)
            for j in range(3):
                p += np.array(self.data.GetPoint(ids.GetId(j)))
            p /= 3.
            p = self.lift_coords(p)
            coords.SetTuple(i, p)
            if i>0:
                cells.InsertNextCell(2)
                cells.InsertCellPoint(i-1)
                cells.InsertCellPoint(i)
        pts = vtk.vtkPoints()
        pts.SetData(coords)
        self.path_poly.SetPoints(pts)
        self.path_poly.SetLines(cells)
        self.path_actor.VisibilityOn()

    def highlight(self, cellid, poly):
        cell = self.data.GetCell(cellid)
        coords = vtk.vtkFloatArray()
        coords.SetNumberOfComponents(3)
        coords.SetNumberOfTuples(3)
        tri = vtk.vtkCellArray()
        tri.InsertNextCell(3)
        idlist = cell.GetPointIds()
        ids = [ idlist.GetId(0), idlist.GetId(1), idlist.GetId(2) ]
        print(f'ids are {ids}')
        for i, id in enumerate(ids):
            coords.SetTuple(i, self.lift_coords(self.data.GetPoint(id)))
            tri.InsertCellPoint(i)
        pts = vtk.vtkPoints()
        pts.SetData(coords)
        poly.SetPoints(pts)
        poly.SetPolys(tri)
        if self.source != -1 and self.target != -1:
            try:
                path = nx.shortest_path(self.G, source=self.source, target=self.target)
                print(path)
                self.draw_path(path)
            except nx.exception.NetworkXNoPath as e: 
                print(f'no path between cell {self.source} and cell {self.target}')

    def get_cell(self):
        pos = self.GetInteractor().GetEventPosition()
        self.picker.Pick(pos[0], pos[1], 0, self.GetDefaultRenderer())
        return self.picker.GetCellId()

    def right_button_press_event(self, obj, event):
        cellid = self.get_cell()
        print(f'cell is {cellid}')
        if self.source >= 0 and cellid == self.source:
            # unselecting source
            self.source = -1
            self.path_actor.VisibilityOff()
            self.src_actor.VisibilityOff()
        elif self.target >= 0 and cellid == self.target:
            # unselecting target
            self.target = -1
            self.trg_actor.VisibilityOff()
            self.path_actor.VisibilityOff()
        elif self.source == -1:
            # missing source selection
            self.source = cellid 
            self.src_actor.VisibilityOn()
            self.highlight(cellid, self.src_poly)
        else:
            # either missing target selection or both source and target are already defined
            self.target = cellid
            self.trg_actor.VisibilityOn()
            self.highlight(cellid, self.trg_poly)
        self.OnRightButtonDown()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze ridge surface connectivity')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input dataset')
    parser.add_argument('-g', '--graph', type=str, default='graph.pickle', help='Name of file containing connectivity graph')    

    args = parser.parse_args()

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(args.input)
    reader.Update()

    ridge = reader.GetOutput()

    minmax = ridge.GetPointData().GetScalars().GetRange()
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(reader.GetOutputPort())
    ctf = vtk.vtkColorTransferFunction()
    ctf.AddRGBPoint(minmax[0], 1, 1, 0)
    ctf.AddRGBPoint(minmax[1], 0, 0, 1)
    mapper.SetLookupTable(ctf)
    mapper.ScalarVisibilityOn()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    window = vtk.vtkRenderWindow()
    window.AddRenderer(renderer)
    window.SetSize(1920, 1080)
    renderer.ResetCamera()
    interactor = vtk.vtkRenderWindowInteractor()
    style = MouseInteractorStyle(ridge, args.graph, renderer)
    style.SetDefaultRenderer(renderer)
    interactor.SetInteractorStyle(style)
    interactor.SetRenderWindow(window)
    interactor.Initialize()
    window.Render()
    interactor.Start()