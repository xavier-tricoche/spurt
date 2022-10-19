#ifndef __VTK_IMAGE_HELPER_HPP__
#define __VTK_IMAGE_HELPER_HPP__

#include "vtk_macros.hpp"
#include "vtk_data_helper.hpp"

#include <third_party/fastmathparser/exprtk.hpp>

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
#include <vtkNrrdReader.h>
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

#if __VTK_HAS_LIC__
// #include <vtkLineIntegralConvolution2D.h>
#endif

namespace vtk_utils {

inline vtkOrientationMarkerWidget* coord_axes_widget(
    vtkRenderWindowInteractor* inter,
    const nvis::vec3& txt_col=nvis::vec3(1,1,1),
    const nvis::vec3& bg_color=nvis::vec3(0,0,0),
    const nvis::vec4& viewport=nvis::vec4(0,0,0.4,0.4),
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

template<typename T>
inline vtkImageData* expr_to_image(const std::vector<std::string>& expr_str,
                                   const nvis::bbox3& bounds,
                                   const nvis::ivec3& res) {

    typedef T                                              value_t;
    typedef nvis::fixed_vector<value_t, 3>                  vec3_t;
    typedef typename vtk_array_traits<value_t>::array_type array_t;
    typedef exprtk::symbol_table<double >           symbol_table_t;
    typedef exprtk::expression<double>                expression_t;
    typedef exprtk::parser<double>                        parser_t;


    const size_t deg=expr_str.size();

    symbol_table_t symbol_table;
    std::vector<expression_t> expr(deg);
    double x, y, z;

    symbol_table.add_variable("x", x);
    symbol_table.add_variable("y", y);
    symbol_table.add_variable("z", z);
    symbol_table.add_constants();

    // std::cout << "symbol table created\n";

    for (int i=0; i<deg; ++i) {
        expr[i].register_symbol_table(symbol_table);
    }

    // std::cout << "symbol table registered to " << deg << " expressions\n";

    parser_t parser;
    for (int i=0; i<deg; ++i) {
        if (!parser.compile(expr_str[i], expr[i])) {
            std::cerr << "parser error: " << parser.error()
                << " for expression " << expr_str[i] << '\n';
            throw std::runtime_error("Invalid expression: " + expr_str[i]);
        }
    }

    // std::cout << "res=" << res << '\n';

    vtkSmartPointer<array_t> values(array_t::New());
    values->SetNumberOfComponents(deg);
    if (deg==1) values->SetNumberOfTuples(res[0]*res[1]*res[2]);
    nvis::vec3 delta=bounds.size()/(nvis::vec3(res)-nvis::vec3(1.));

    if (res[2]==1) delta[2]=1.;

    // std::cout << "delta=" << delta << '\n';

    int counter=0;
    for (int n=0; n<res[0]*res[1]*res[2]; ++n) {
        // n = i + res[0]*(j + res[1]*k)
        int i=n%res[0];
        int j=(n/res[0])%res[1];
        int k=n/(res[0]*res[1]);
        x=std::min(bounds.min()[0]+i*delta[0], bounds.max()[0]);
        y=std::min(bounds.min()[1]+j*delta[1], bounds.max()[1]);
        z=std::min(bounds.min()[2]+k*delta[2], bounds.max()[2]);
        // std::cout << "evaluation at " << nvis::vec3(x,y,z) << '\n';
        if (deg==1) {
            value_t v=static_cast<value_t>(expr[0].value());
            values->SetValue(counter, v);
        }
        else if (deg==3) {
            vec3_t v;
            v[0]=static_cast<value_t>(expr[0].value());
            v[1]=static_cast<value_t>(expr[1].value());
            v[2]=static_cast<value_t>(expr[2].value());
            // std::cout << "v=" << v << '\n';
            values->InsertNextTypedTuple(v.begin());
        }
    }

    /*
    for (z=bounds.min()[2]; z<=bounds.max()[2]; z+=delta[2]) {
        for (y=bounds.min()[1]; y<=bounds.max()[1]; y+=delta[1]) {
            for (x=bounds.min()[0]; x<=bounds.max()[0]; x+=delta[0], ++counter) {
                std::cout << "evaluation at " << nvis::vec3(x,y,z) << '\n';
                if (deg==1) {
                    double v=expr[0].value();
                    values->SetValue(counter, v);
                }
                else if (deg==3) {
                    nvis::vec3 v;
                    v[0]=expr[0].value();
                    v[1]=expr[1].value();
                    v[2]=expr[2].value();
                    std::cout << "v=" << v << '\n';
                    values->InsertNextTypedTuple(v.begin());
                }
            }
        }
    }
    */

    vtkImageData* out=vtkImageData::New();
    out->SetDimensions(res[0], res[1], res[2]);
    out->SetOrigin(bounds.min()[0], bounds.min()[1], bounds.min()[2]);
    out->SetSpacing(delta[0], delta[1], delta[2]);
    if (deg==1) {
        out->GetPointData()->SetScalars(values);
    }
    else if (deg==3) {
        out->GetPointData()->SetVectors(values);
    }
    return out;
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
        writer->SetInputData(tmp_img);
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

template<>
inline void save_image_in_given_format<vtkTIFFWriter>(const vtkImageData* data,
                                                      const std::string& filename,
													  int quality) {
    VTK_CREATE(vtkTIFFWriter, writer);
    VTK_CONNECT(writer, const_cast<vtkImageData*>(data));
    writer->SetFileName(filename.c_str());
    if (quality < 100) {
        writer->SetCompressionToJPEG();
    }
    else {
        writer->SetCompressionToPackBits();
    }
    writer->Write();
}

inline vtkImageData* load_image(const std::string& filename)
{
    std::string ext = xavier::filename::extension(filename);
    if (ext == "png") {
        return load_image_of_given_format<vtkPNGReader>(filename);
    } else if (ext == "tif" || ext == "tiff") {
        return load_image_of_given_format<vtkTIFFReader>(filename);
    } else if (ext == "bmp") {
        return load_image_of_given_format<vtkBMPReader>(filename);
    } else if (ext == "jpg" || ext == "jpeg" || ext == "jpg2" || ext == "jpeg2") {
        return load_image_of_given_format<vtkJPEGReader>(filename);
    }
    else if (ext == "nrrd" || ext == "nhdr") {
        return load_image_of_given_format<vtkNrrdReader>(filename);
    }
    else if (ext == "pnm" || ext == "ppm") {
        return load_image_of_given_format<vtkPNMReader>(filename);
    } else {
        std::cerr << "image type of " << filename << " (" << ext << ") unrecognized\n";
        throw std::runtime_error("unrecognized image filename extension");
    }
}

inline void save_image(const vtkImageData* data, const std::string& filename, int quality=100)
{
    std::string ext = xavier::filename::extension(filename);
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
        save_image_in_given_format<vtkTIFFWriter>(data, filename, quality);
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

#if 0 /*__VTK_HAS_LIC__*/
// inline vtkImageData* do_lic(const vtkImageData* rhs, int nsteps, float eps,
//                             bool enhanced=false, int factor=1)
// {
//     // temporary render window provides local context to GPU computation
//     VTK_PTR(vtkRenderWindow, context);
//
//     // compute LIC texture
//     vtkImageData* img;
//     VTK_PTR(vtkLineIntegralConvolution2D, lic);
//     VTK_CONNECT(lic, const_cast<vtkImageData*>(rhs));
//     lic->SetContext(context);
//     lic->SetNumberOfSteps(nsteps);
//     lic->SetStepSize(eps);
//     lic->SetMagnification(factor);
// 	lic->EnhancedLICOn();
// 	lic->AntiAliasOn();
//     // if (enhanced) lic->EnhanceOn();
//     lic->UpdateInformation();
//     lic->Update();
//     context->Delete();
//     img=lic->GetOutput();
//     img->Register(img);
//     lic->Delete();
//     return img;
// }
#endif

inline vtkRenderer* fill_window(vtkRenderer* inout, const nvis::bbox2& bounds) {
    inout->GetActiveCamera()->SetParallelProjection(1);
    nvis::vec2 center=bounds.center();
    inout->GetActiveCamera()->SetPosition(center[0], center[1], 1);
    inout->GetActiveCamera()->SetFocalPoint(center[0], center[1], 0);
    inout->GetActiveCamera()->SetViewUp(0, 1, 0);
    inout->GetActiveCamera()->SetParallelScale(0.5*(bounds.max()[1]-bounds.min()[1]));
    return inout;
}

} // vtk_utils

#endif
