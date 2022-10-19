#ifndef __VTK_COLORBAR_HELPER_HPP__
#define __VTK_COLORBAR_HELPER_HPP__

#include <string>
#include <fstream>
#include <math/fixed_vector.hpp>

#include <vtkColorTransferFunction.h>
#include <vtkProperty2D.h>
#include <vtkScalarBarActor.h>
#include <vtkSmartPointer.h>
#include <vtkTextProperty.h>

namespace vtk_utils {
    
typedef spurt::vec3 color_t;

void export_colormap(const std::string& name, vtkColorTransferFunction* ctf) {
    std::fstream out(name.c_str(), std::ios::out);
    int npts = ctf->GetSize();
    double val[6];
    for (int i=0; i<npts; ++i) {
        ctf->GetNodeValue(i, val);
        out << val[0] << " " << val[1] << " " << val[2] << " " << val[3] << '\n';
    }
    out.close();
}

void import_colormap(const std::string& name, vtkColorTransferFunction* ctf) {
    std::fstream in(name.c_str(), std::ios::in);
    double v, r, g, b;
    while (!in.eof()) {
        in >> v >> r >> g >> b;
        if (in.good()) ctf->AddRGBPoint(v, r, g, b);
    }
    in.close();
}
    
struct colorbar_param {
    colorbar_param() 
        : title("No Title"), title_col(1.,1.,1.), label_col(1.,1.,1.),
          pos(0.9, 0.5), 
          width(80), height(400),
          nlabels(4), font_size(22),
          title_offset(20) {}
  
    std::string   title;
    color_t       title_col, label_col;
    spurt::vec2    pos;
    unsigned int  width, height;
    unsigned int  nlabels;
    unsigned int  font_size;
    int           title_offset;
};

vtkSmartPointer<vtkScalarBarActor> 
colorbar(vtkSmartPointer<vtkColorTransferFunction> ctf,
         const colorbar_param& param = colorbar_param(), 
         bool is_float=true)
{
    // Create a color bar
    vtkSmartPointer<vtkScalarBarActor> scalar_bar =
        vtkSmartPointer<vtkScalarBarActor>::New();
    
    // size and relative position
    scalar_bar->SetLookupTable(ctf);
    scalar_bar->SetPosition(param.pos[0], param.pos[1]);
    scalar_bar->SetMaximumWidthInPixels(param.width);
    scalar_bar->SetMaximumHeightInPixels(param.height);
    
    // title properties
    scalar_bar->SetTitle(param.title.c_str());
    scalar_bar->GetTitleTextProperty()->SetColor(param.title_col[0], 
                                                 param.title_col[1], 
                                                 param.title_col[2]);
    scalar_bar->GetTitleTextProperty()->SetLineOffset(param.title_offset);
    scalar_bar->GetTitleTextProperty()->ShadowOff();
    scalar_bar->GetTitleTextProperty()->SetFontSize(param.font_size);
    scalar_bar->GetTitleTextProperty()->BoldOn();
    
    // label properties
    scalar_bar->SetNumberOfLabels(param.nlabels);
    scalar_bar->SetLabelFormat(is_float ? "%0.2f" : "%0.0f");
    scalar_bar->GetLabelTextProperty()->SetColor(param.label_col[0], 
                                                 param.label_col[1], 
                                                 param.label_col[2]);
    scalar_bar->GetLabelTextProperty()->SetFontSize(param.font_size);
    scalar_bar->GetLabelTextProperty()->BoldOff();
    scalar_bar->GetLabelTextProperty()->ShadowOff();
    
    return scalar_bar;
}

} // vtk_utils

#endif
