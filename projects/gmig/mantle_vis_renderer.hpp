/*=========================================================================

  Copyright 2004 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/

/*========================================================================
 For general information about using VTK and Qt, see:
 http://www.trolltech.com/products/3rdparty/vtksupport.html
=========================================================================*/

/*========================================================================
 !!! WARNING for those who want to contribute code to this file.
 !!! If you use a commercial edition of Qt, you can modify this code.
 !!! If you use an open source version of Qt, you are free to modify
 !!! and use this code within the guidelines of the GPL license.
 !!! Unfortunately, you cannot contribute the changes back into this
 !!! file.  Doing so creates a conflict between the GPL and BSD-like VTK
 !!! license.
=========================================================================*/

#ifndef _ALT_MANTLE_VIS_HPP_
#define _ALT_MANTLE_VIS_HPP_

#include <QMainWindow>
#include "ui_mantle_vis_renderer.h"
#include <vtkSmartPointer.h>

class vtkActor;
class vtkRenderWindow;
class vtkObject;
class vtkDataSetMapper;
class vtkOrientationMarkerWidget;
class vtkColorTransferFunction;
class vtkPolyData;

class QAction;

// Forward Qt class declarations
class Ui_MainWindow;

namespace spurt { namespace gmig {
    
struct TypeAttributes;

class mantle_vis_control;

class mantle_vis_renderer : public QMainWindow
{
    Q_OBJECT
public:
    mantle_vis_renderer(int argc, char* argv[]);
    ~mantle_vis_renderer();

private slots:
    // slot functions for GUI signals
    void slot_file_action(QAction*);
    
    // slot functions for signals emitted by control window UI
    void slot_points_check(int);
    void slot_edges_check(int);
    void slot_cells_check(int);
    void slot_axes_check(int);
    void slot_world_check(int);
    void slot_plus_check(int);
    void slot_minus_check(int);
    void slot_radius_slider(int);
    void slot_radius_spinbox(double);
    void slot_size_spinbox(int);
    void slot_min_spinbox(int);
  
protected:    
    void load_file(const std::string&);
    void save_snapshot();
    void draw();
    void update(bool=true);
    void internal_main();
    void camera_out();
    void camera_in();
    void show_hide_control(QAction*);
    void save_patches(TypeAttributes&);
    
    void hide_points(TypeAttributes&, bool=true);
    void hide_points(bool=true);
    void show_points(TypeAttributes&, bool=true);
    void show_points(bool=true);
    void hide_edges(TypeAttributes&, bool=true);
    void hide_edges(bool=true);
    void show_edges(TypeAttributes&, bool=true);
    void show_edges(bool=true);
    void hide_cells(TypeAttributes&, bool=true);
    void hide_cells(bool=true);
    void show_cells(TypeAttributes&, bool=true);
    void show_cells(bool=true);
    void redraw_points(TypeAttributes&, bool=true);
    void redraw_edges(TypeAttributes&, bool=true);
    void redraw_cells(TypeAttributes&, bool=true);
    bool slot_type_check(int, TypeAttributes&);
    void draw_from_scratch(TypeAttributes&, bool=true, bool=false);
    
    mantle_vis_control *controls;
    
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkOrientationMarkerWidget> axes_widget;
    vtkSmartPointer<vtkColorTransferFunction> color_tf;
    
    vtkSmartPointer<vtkPolyData>
        assign_component_ids(vtkSmartPointer<vtkPolyData>, const std::vector< std::vector<int> >&, bool=false);
    vtkSmartPointer<vtkPolyData>
        assign_value_by_component(vtkSmartPointer<vtkPolyData>, const std::vector< std::vector<int> >&);
    
private:

    // Designer form
    Ui_MainWindow *ui;
    bool initialized;
};

} // gmig
} // spurt

#endif // _ALT_MANTLE_VIS_HPP_


