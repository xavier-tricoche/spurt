#ifndef __VTK_IO_HELPER_HPP__
#define __VTK_IO_HELPER_HPP__

#include <fstream>
#include <string>
#include <vector>

#include "vtk_macros.hpp"
#include <format/filename.hpp>
#include <misc/strings.hpp>

#include <math/fixed_vector.hpp>

#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkRectilinearGrid.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkSliderRepresentation.h>
#include <vtkSliderRepresentation2D.h>

#include <vtkXMLImageDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
    

namespace vtk_utils {
    
template< typename Reader_ >
vtkDataSet* __readVTK(const std::string& name) {
    typedef Reader_ reader_type;
    
    vtkSmartPointer<reader_type> reader = vtkSmartPointer<reader_type>::New();
    reader->SetFileName(name.c_str());
    reader->Update();
    vtkDataSet* data = reader->GetOutput();
    data->Register(0);
    return data;
}

template< typename DataSetPtr_, typename Writer_ >
void __saveVTK(DataSetPtr_ data, const std::string& name) {
    typedef Writer_ writer_type;
    
    // note: no check is performed that data type is compatible with selected writer type
    vtkSmartPointer<writer_type> writer = vtkSmartPointer<writer_type>::New();
    writer->SetFileName(name.c_str());
    writer->SetInputData(data);
    writer->Write();
} 

vtkDataSet* readVTK(const std::string& filename)
{
    std::string ext = spurt::filename::extension(filename);
    spurt::lower_case(ext);
    if      (ext == "vtk") return __readVTK<vtkDataSetReader>(filename);
    else if (ext == "vti") return __readVTK<vtkXMLImageDataReader>(filename);
    else if (ext == "vtu") return __readVTK<vtkXMLUnstructuredGridReader>(filename);
    else if (ext == "vtp") return __readVTK<vtkXMLPolyDataReader>(filename);
    else throw std::runtime_error("Unrecognized VTK file extension");
    
    return NULL;
}

template< typename DataSetPtr_ >
void saveVTK(DataSetPtr_ dataset, const std::string& filename) {
    std::string ext = spurt::filename::extension(filename);
    spurt::lower_case(ext);
    if      (ext == "vtk") __saveVTK<DataSetPtr_, vtkDataSetWriter>(dataset, filename);
    else if (ext == "vti") __saveVTK<DataSetPtr_, vtkXMLImageDataWriter>(dataset, filename);
    else if (ext == "vtu") __saveVTK<DataSetPtr_, vtkXMLUnstructuredGridWriter>(dataset, filename);
    else if (ext == "vtp") __saveVTK<DataSetPtr_, vtkXMLPolyDataWriter>(dataset, filename);
    else if (ext == "vts") __saveVTK<DataSetPtr_, vtkXMLStructuredGridWriter>(dataset, filename);
    else if (ext == "vtr") __saveVTK<DataSetPtr_, vtkXMLRectilinearGridWriter>(dataset, filename);
    else throw std::runtime_error("Unrecognized VTK file extension");
}

template< typename DataSetPtr_ >
void saveVTK_XML(DataSetPtr_ dataset, const std::string& filename) {
    if (vtkImageData::SafeDownCast(dataset)) {
        spurt::filename::replace_extension(filename, "vti");
    }
    else if (vtkUnstructuredGrid::SafeDownCast(dataset)) {
        spurt::filename::replace_extension(filename, "vtu");
    }
    else if (vtkPolyData::SafeDownCast(dataset)) {
        spurt::filename::replace_extension(filename, "vtp");
    }
    else if (vtkRectilinearGrid::SafeDownCast(dataset)) {
        spurt::filename::replace_extension(filename, "vtr");
    }
    else if (vtkStructuredGrid::SafeDownCast(dataset)) {
        spurt::filename::replace_extension(filename, "vts");
    }
    else {
        spurt::filename::replace_extension(filename, "vtk");
        std::cout << "WARNING: Unrecognized VTK dataset type. Using Legacy format" << '\n';
    }
    
    saveVTK(dataset, filename);
};

// vtkDataSet* readDLR(const std::string& gridname, const std::string& pvalname)
// {
//     vtkSmartPointer<vtkDLRReader> reader = vtkSmartPointer<vtkDLRReader>::New();
//     reader->SetGridFileName(grid.c_str());
//     reader->SetPvalFileName(pval.c_str());
//     reader->Update();
//     
//     vtkDataSet* grid = static_cast<vtkDataSet*>(reader->GetOutput());
//     grid->Register(0);
//     reader->Delete();
//     return grid;
// }

void save_lines_VTK(const std::string& name, const std::string& comment,
                    const std::vector<std::vector<spurt::vec3> >& lines)
{
    int npts = 0;
    for (int i = 0 ; i < lines.size() ; ++i) {
        npts += lines[i].size();
    }
    
    std::fstream file(name.c_str(), std::ios::out);
    file << "# vtk DataFile Version 2.0\n"
         << comment << '\n'
         << "ASCII\n"
         << "DATASET POLYDATA\n"
         << "POINTS " << npts  << " float\n";
    for (int i = 0 ; i < lines.size() ; ++i) {
        for (int j = 0 ; j < lines[i].size() ; ++j) {
            const spurt::vec3& x = lines[i][j];
            file << x[0] << " " << x[1] << " " << x[2] << '\n';
        }
    }
    file << "LINES " << lines.size() << " " << npts + lines.size() << '\n';
    int c = 0;
    for (int i = 0 ; i < lines.size() ; ++i) {
        file << lines[i].size();
        for (int j = 0 ; j < lines[i].size() ; ++j) {
            file << " " << c++;
        }
        file << '\n';
    }
    
    file.close();
}

void save_triangles_VTK(const std::string& name, const std::string& comment,
                        const std::vector<spurt::vec3>& points,
                        const std::vector<spurt::vec3>& triangles)
{
    int npts = points.size();
    
    std::fstream file(name.c_str(), std::ios::out);
    file << "# vtk DataFile Version 2.0\n"
         << comment << '\n'
         << "ASCII\n"
         << "DATASET POLYDATA\n"
         << "POINTS " << npts  << " float\n";
    for (int i = 0 ; i < points.size() ; ++i) {
        const spurt::vec3& x = points[i];
        file << x[0] << " " << x[1] << " " << x[2] << '\n';
    }
    file << "POLYGONS " << triangles.size() << " " << 4*triangles.size() << '\n';
    for (int i = 0 ; i < triangles.size() ; ++i) {
        const spurt::vec3& t = triangles[i];
        file << "3 " << t[0] << " " << t[1] << " " << t[2] << '\n';
    }
    
    file.close();
}

vtkStructuredPoints* load_nrrd(const std::string& filename) {
    Nrrd* nin = spurt::readNrrd(filename);
    int attribute_type = 0; // scalar
    if (nin->axis[0].size == 2 || nin->axis[0].size == 3) attribute_type = 1;
    else if (nin->axis[0].size == 4 || nin->axis[0].size == 9) attribute_type = 2;
    
    int dims[3] = {1, 1, 1};
    double spc[3] = {1, 1, 1};
    double min[3] = {0, 0, 0};
    int shift = (attribute_type ? 1 : 0);
    for (int i=0 ; i<nin->dim ; ++i) {
        dims[i] = nin->axis[shift+i].size;
        double sp = nin->axis[shift+i].spacing;
        double m = nin->axis[i+shift].min;
        if (!std::isinf(sp) && !std::isnan(sp)) spc[i] = sp;
        if (!std::isinf(m) && !std::isnan(m)) min[i] = m;
    }
    
    VTK_PTR(vtkStructuredPoints, dataset);
    dataset->SetDimensions(dims[0], dims[1], dims[2]);
    dataset->SetSpacing(spc[0], spc[1], spc[2]);
    dataset->SetOrigin(min[0], min[1], min[2]);
    dataset->SetExtent(0, dims[0]-1, 0, dims[1]-1, 0, 0);
    if (nin->type == nrrdTypeFloat) {
        std::vector<float> data;
        spurt::to_vector(data, nin);
        if (attribute_type == 0) 
            add_scalars(dataset, data);
        else if (attribute_type == 1) {
            if (nin->axis[0].size == 2)
                add_vectors_from_numbers<std::vector<float>, vtkStructuredPoints*, 2>(dataset, data);
            else
                add_vectors_from_numbers<std::vector<float>, vtkStructuredPoints*, 3>(dataset, data);
        }
        else if (attribute_type == 2) {
            if (nin->axis[0].size == 4)
                add_tensors_from_numbers<std::vector<float>, vtkStructuredPoints*, 4>(dataset, data);
            else
                add_tensors_from_numbers<std::vector<float>, vtkStructuredPoints*, 9>(dataset, data);
        }
    }
    else if (nin->type == nrrdTypeDouble) {
        std::vector<double> data;
        spurt::to_vector(data, nin);
        if (attribute_type == 0) 
            add_scalars(dataset, data);
        else if (attribute_type == 1) {
            if (nin->axis[0].size == 2)
                add_vectors_from_numbers<std::vector<double>, vtkStructuredPoints*, 2>(dataset, data);
            else
                add_vectors_from_numbers<std::vector<double>, vtkStructuredPoints*, 3>(dataset, data);
        }
        else if (attribute_type == 2) {
            if (nin->axis[0].size == 4)
                add_tensors_from_numbers<std::vector<double>, vtkStructuredPoints*, 4>(dataset, data);
            else
                add_tensors_from_numbers<std::vector<double>, vtkStructuredPoints*, 9>(dataset, data);
        }
    }
    else {
        std::vector<int> data;
        spurt::to_vector(data, nin);
        if (attribute_type == 0) 
            add_scalars(dataset, data);
        else if (attribute_type == 1) {
            if (nin->axis[0].size == 2)
                add_vectors_from_numbers<std::vector<int>, vtkStructuredPoints*, 2>(dataset, data);
            else
                add_vectors_from_numbers<std::vector<int>, vtkStructuredPoints*, 3>(dataset, data);
        }
        else if (attribute_type == 2) {
            if (nin->axis[0].size == 4)
                add_tensors_from_numbers<std::vector<int>, vtkStructuredPoints*, 4>(dataset, data);
            else
                add_tensors_from_numbers<std::vector<int>, vtkStructuredPoints*, 9>(dataset, data);
        }
    }
    nrrdNuke(nin);
    
    return dataset;
}

vtkSliderRepresentation* 
make_slider_representation(const std::string& text, const spurt::vec2& range, 
                           double value, double x, double dx, double y=0.07, 
                           const std::string& format="") {
    VTK_PTR(vtkSliderRepresentation2D, rep);
    rep->SetMinimumValue(range[0]);
    rep->SetMaximumValue(range[1]);
    rep->SetValue(value);
    rep->SetTitleText(text.c_str());
    rep->SetLabelFormat(format.empty() ? "%0.0f" : format.c_str());
    rep->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
    rep->GetPoint1Coordinate()->SetValue(x, y);
    rep->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
    rep->GetPoint2Coordinate()->SetValue(x+dx, y);
    rep->SetSliderLength(0.005);
    rep->SetSliderWidth(0.02);
    rep->SetEndCapLength(0.002);
    rep->SetEndCapWidth(0.03);
    rep->SetTubeWidth(0.002);
    rep->GetTubeProperty()->SetColor(0.9, 0.9, 0.9);
    rep->GetCapProperty()->SetColor(0.9, 0.9, 0.9);
    rep->GetSliderProperty()->SetColor(1, 0, 0);
    rep->GetTitleProperty()->SetColor(1, 1, 1);
    rep->GetTitleProperty()->ShadowOff();
    rep->GetTitleProperty()->SetFontFamilyToTimes();
    rep->SetTitleHeight(0.02);
    rep->GetTitleProperty()->BoldOff();
    rep->GetLabelProperty()->SetColor(1, 1, 1);
    rep->SetLabelHeight(0.02);
    rep->GetLabelProperty()->SetFontFamilyToTimes();
    rep->GetLabelProperty()->BoldOff();
    rep->GetLabelProperty()->ShadowOff();
    return rep;
}

} // vtk_utils


#endif



