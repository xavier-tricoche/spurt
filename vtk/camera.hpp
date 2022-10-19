#ifndef __CS530_CAMERA_HPP__
#define __CS530_CAMERA_HPP__

// STL
#include <iostream>
#include <string>
// Qt
#include <QMenu>
#include <QFileDialog>
#include <QString>
#include <QWidget>
// VTK
#include "vtkCamera.h"
#include "vtkRenderer.h"
// nvis
#include <math/fixed_vector.hpp>

namespace vtk_utils {

inline
void export_camera_settings(QWidget* parent, vtkRenderer* renderer,
                            const std::string& base_dir) {
    QString fn = 
        QFileDialog::getSaveFileName(parent, 
                                     QString("Select an output camera file"),
                                     QString(base_dir.c_str()));
    std::string filename = fn.toStdString();
    export_camera_settings(filename, renderer);
};

inline
void import_camera_settings(QWidget* parent, vtkRenderer* renderer,
                            const std::string& base_dir) {
    QString fn = 
        QFileDialog::getOpenFileName(parent,
                                     QString("Select an input camera file"),
                                     QString(base_dir.c_str()));
    std::string filename = fn.toStdString();
    import_camera_settings(filename, renderer);

} // vtk_utils

#endif
