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

#ifndef _MANTLE_VIS_CONTROL_HPP_
#define _MANTLE_VIS_CONTROL_HPP_

#include <QMainWindow>
#include <ui_mantle_vis_control.h>
#include <vtkSmartPointer.h>

class vtkActor;
class vtkRenderWindow;
class vtkObject;
class vtkDataSetMapper;

class QAction;

// Forward Qt class declarations
class Ui_ControlWidget;

namespace xavier { namespace gmig {

class mantle_vis_control : public QWidget
{    
    Q_OBJECT
public:
    mantle_vis_control();
    ~mantle_vis_control();

private slots:
    // slot functions for Qt signals
    // void slot_points_check(int);
    // void slot_edges_check(int);
    // void slot_axes_check(int);
    // void slot_world_check(int);
    // void slot_plus_check(int);
    // void slot_minus_check(int);
    // void slot_radius_slider(double);
    // void slot_radius_spinbox(double);
    // void slot_size_spinbox(int);
    
signals:
    // signals for communication with parent render window
    // void points_status_changed(int on);
    // void edges_status_changed(int on);
    // void axes_status_changed(int on);
    // void world_status_changed(int on);
    // void plus_status_changed(int on);
    // void minus_status_changed(int on);
    // void neigh_radius_changed(double r);
    // void neigh_size_changed(int n);
      
private:

public:
  // Designer form
  Ui_ControlWidget *ui;
};

} // gmig
} // xavier

#endif // _MANTLE_VIS_CONTROL_HPP_


