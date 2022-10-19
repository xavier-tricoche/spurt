#ifndef __VTK_IMAGE_HELPER_HPP__
#define __VTK_IMAGE_HELPER_HPP__

#include "vtk_macros.hpp"
#include "vtk_data_helper.hpp"
#include <string>

#include <math/fixed_vector.hpp>

#include <vtkAxesActor.h>
#include <vtkBMPReader.h>
#include <vtkBMPWriter.h>
#include <vtkCaptionActor2D.h>
#include <vtkDataArray.h>
#include <vtkDataSetWriter.h>
#include <vtkImageData.h>
#include <vtkImageShiftScale.h>
#include <vtkJPEGReader.h>
#include <vtkJPEGWriter.h>
#if VTK_MAJOR_VERSION >= 6
#include <vtkNrrdReader.h>
#endif
#include <vtkOrientationMarkerWidget.h>
#include <vtkPNGReader.h>
#include <vtkPNGWriter.h>
#include <vtkPNMReader.h>
#include <vtkPNMWriter.h>
#include <vtkPointData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTexture.h>
#include <vtkTIFFReader.h>
#include <vtkTIFFWriter.h>
#include <vtkWindowToImageFilter.h>

namespace vtk_utils {

inline vtkOrientationMarkerWidget* coord_axes_widget(
    vtkRenderWindowInteractor* inter,
    const spurt::vec3& txt_col=spurt::vec3(1,1,1),
    const spurt::vec3& bg_color=spurt::vec3(0,0,0),
    const spurt::vec4& viewport=spurt::vec4(0,0,0.4,0.4),
    double radius=0.05,
    const std::string& Xlabel="X",
    const std::string& Ylabel="Y",
    const std::string Zlabel="Z") {
    VTK_CREATE(vtkAxesActor, axes);
    axes->SetCylinderRadius(radius);
    axes->SetShaftTypeToCylinder();
    axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetColor(txt_col[0], txt_col[1], txt_col[2]);
    axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetColor(txt_col[0], txt_col[1], txt_col[2]);
    axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetColor(txt_col[0], txt_col[1], txt_col[2]);
    axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetBackgroundColor(bg_color[0], bg_color[1], bg_color[2]);
    axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetBackgroundColor(bg_color[0], bg_color[1], bg_color[2]);
    axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->
        SetBackgroundColor(bg_color[0], bg_color[1], bg_color[2]);
    axes->SetXAxisLabelText(Xlabel.c_str());
    axes->SetYAxisLabelText(Ylabel.c_str());
    axes->SetZAxisLabelText(Zlabel.c_str());

    vtkOrientationMarkerWidget* widget=vtkOrientationMarkerWidget::New();
    widget->SetOutlineColor(bg_color[0], bg_color[1], bg_color[2]);
    widget->SetOrientationMarker(axes);
    widget->SetInteractor(inter);
    widget->SetViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
    return widget;
}

inline vtkTexture* image_to_texture(const vtkImageData* img)
{
    const vtkDataArray* data = const_cast<vtkImageData*>(img)->GetPointData()->GetScalars();
    VTK_PTR(vtkTexture, texture);
    if (const_cast<vtkDataArray*>(data)->GetDataType() == VTK_UNSIGNED_CHAR ||
        const_cast<vtkDataArray*>(data)->GetDataType() == VTK_UNSIGNED_SHORT) {
        VTK_IMAGE_CONNECT(texture, const_cast<vtkImageData*>(img));
    } else {
        vtkDataArray* byte_data = convert<unsigned char>(data);
        VTK_CREATE(vtkImageData, tmp_img);
        tmp_img->CopyStructure(const_cast<vtkImageData*>(img));
        tmp_img->GetPointData()->SetScalars(byte_data);
        //
        VTK_PTR(vtkDataSetWriter, writer);
#if (VTK_MAJOR_VERSION >= 6)
        writer->SetInputData(tmp_img);
#else
        writer->SetInput(tmp_img);
#endif
        writer->SetFileName("debug2.vtk");
        writer->SetFileTypeToBinary();
        writer->Update();
        //
        byte_data->Delete();
        VTK_IMAGE_CONNECT(texture, tmp_img);
    }
    texture->InterpolateOn();
    return texture;
}

template<typename Reader>
inline vtkImageData* load_image_of_given_format(const std::string& filename)
{
    VTK_CREATE(Reader, reader);
    reader->SetFileName(filename.c_str());
    reader->Update();
    vtkImageData* _img = vtkImageData::SafeDownCast(reader->GetOutput());
    reader->GetOutput()->Register(_img);
    return _img;
    reader->Delete();
}

template<typename Writer>
inline void save_image_in_given_format(const vtkImageData* data,
                                       const std::string& filename,
                                       int quality=100)
{
    VTK_CREATE(Writer, writer);
    VTK_CONNECT(writer, const_cast<vtkImageData*>(data));
    writer->SetFileName(filename.c_str());
    writer->Write();
}

template<>
inline void save_image_in_given_format<vtkJPEGWriter>(const vtkImageData* data,
                                                      const std::string& filename,
                                                      int quality) {
    VTK_CREATE(vtkJPEGWriter, writer);
    VTK_CONNECT(writer, const_cast<vtkImageData*>(data));
    writer->SetFileName(filename.c_str());
    writer->SetQuality(quality);
    writer->Write();
}

inline vtkImageData* load_image(const std::string& filename)
{
    std::string ext = spurt::filename::extension(filename);
    if (ext == "png") {
        return load_image_of_given_format<vtkPNGReader>(filename);
    } else if (ext == "tif" || ext == "tiff") {
        return load_image_of_given_format<vtkTIFFReader>(filename);
    } else if (ext == "bmp") {
        return load_image_of_given_format<vtkBMPReader>(filename);
    } else if (ext == "jpg" || ext == "jpeg" || ext == "jpg2" || ext == "jpeg2") {
        return load_image_of_given_format<vtkJPEGReader>(filename);
    }
#if VTK_MAJOR_VERSION >= 6
    else if (ext == "nrrd" || ext == "nhdr") {
        return load_image_of_given_format<vtkNrrdReader>(filename);
    }
#endif
    else if (ext == "pnm" || ext == "ppm") {
        return load_image_of_given_format<vtkPNMReader>(filename);
    } else {
        std::cerr << "image type of " << filename << " (" << ext << ") unrecognized\n";
        throw std::runtime_error("unrecognized image filename extension");
    }
}

inline void save_image(const vtkImageData* data, const std::string& filename, int quality=100)
{
    std::string ext = spurt::filename::extension(filename);
    if (ext == "png") {
        VTK_CREATE(vtkImageShiftScale, cast);
        VTK_CONNECT(cast, const_cast<vtkImageData*>(data));
        double* range = const_cast<vtkImageData*>(data)->GetPointData()->GetScalars()->GetRange();
        double _min = range[0];
        double _max = range[1];
        cast->SetShift(-_min);
        cast->SetScale(255./(_max-_min));
        cast->SetOutputScalarTypeToUnsignedChar();
        cast->Update();
        save_image_in_given_format<vtkPNGWriter>(cast->GetOutput(), filename);
    } else if (ext == "tif" || ext == "tiff") {
        save_image_in_given_format<vtkTIFFWriter>(data, filename);
    } else if (ext == "bmp") {
        save_image_in_given_format<vtkBMPWriter>(data, filename);
    } else if (ext == "jpg" || ext == "jpeg" || ext == "jpg2" || ext == "jpeg2") {
        save_image_in_given_format<vtkJPEGWriter>(data, filename, quality);
    }
#if __VTK_NRRD_WRITER_WORKS__
    else if (ext == "nrrd" || ext == "nhdr") {
        save_image_in_given_format<vtkNrrdWriter>(data, filename);
    }
#endif
    else if (ext == "pnm" || ext == "ppm") {
        save_image_in_given_format<vtkPNMWriter>(data, filename);
    } else if (ext == "vtk") {
        save_image_in_given_format<vtkDataSetWriter>(data, filename);
    } else if (ext.empty()) {
        save_image_in_given_format<vtkPNGWriter>(data, filename + ".png");
    } else {
        std::cerr << "image type of " << filename << " (" << ext << ") unrecognized\n";
        throw std::runtime_error("unrecognized image filename extension");
    }
}

inline void save_frame(vtkRenderWindow* window,
                       const std::string& filename,
                       int quality=100)
{
    // turn screen (frame buffer) content into an image
    VTK_CREATE(vtkWindowToImageFilter, capture);
    capture->SetInput(window);
    capture->Update();
    save_image(capture->GetOutput(), filename, quality);
}


inline vtkRenderer* fill_window(vtkRenderer* inout, const spurt::bbox2& bounds) {
    inout->GetActiveCamera()->SetParallelProjection(1);
    spurt::vec2 center=bounds.center();
    inout->GetActiveCamera()->SetPosition(center[0], center[1], 1);
    inout->GetActiveCamera()->SetFocalPoint(center[0], center[1], 0);
    inout->GetActiveCamera()->SetViewUp(0, 1, 0);
    inout->GetActiveCamera()->SetParallelScale(0.5*(bounds.max()[1]-bounds.min()[1]));
    return inout;
}

} // vtk_utils

#endif
