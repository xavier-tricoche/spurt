#include <iterator>
#include <map>
#include <memory>
#include <string>

// xavier's utilities
#include <format/format.hpp>
#include <misc/option_parse.hpp>

// VTK
#include <vtkCaptionActor2D.h>
#include <vtkColorTransferFunction.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

// VTK helper functions and macros
#include <VTK/vtk_utils.hpp>
#include <graphics/colors.hpp>

// Qt
#include <QMenu>
#include <QFileDialog>
#include "QVTKWidget.h"
#include "QVTKInteractor.h"
#include "main_window.hpp"
#include "mantle_vis_renderer.hpp"

// Tapkee
// #include <tapkee/tapkee.hpp>
// #include <tapkee/methods.hpp>

// miscellaneous
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

std::string home_directory;
vtkSmartPointer< vtkActor > iso0_actor, iso1_actor;


// Qt slots
void corn_mri_main_window::slot_set_value0(double v) {

}
void corn_mri_main_window::slot_set_value1(double v) {

}
void corn_mri_main_window::slot_set_min(double v) {

}
void corn_mri_main_window::slot_set_max(double v) {

}
void corn_mri_main_window::slot_show_value0(int state) {
	if (state == Qt::Checked) {
		show_isosurface(0);
	}
	else {
		hise_isosurface(0);
	}
}
void corn_mri_main_window::slot_show_value1(int state) {
	if (state == Qt::Checked) {
		show_isosurface(1);
	}
	else {
		hise_isosurface(1);
	}
}
void corn_mri_main_window::slot_set_label(Qstring txt) {

}
void corn_mri_main_window::slot_set_id(QString txt) {

}
void corn_mri_main_window::slot_front_view() {

}
void corn_mri_main_window::slot_back_view() {

}
void corn_mri_main_window::slot_top_push() {

}
void corn_mri_main_window::slot_bottom_view() {

}
void corn_mri_main_window::slot_left_view() {

}
void corn_mri_main_window::slot_right_view() {

}
void corn_mri_main_window::slot_compute() {

}

corn_mri_main_window::corn_mri_main_window(int argc, char* argv[]) {

	// initialization
	this->ui = new main_window;
	this->ui->setupUi(this);

	VTK_CREATE(Renderer, main_renderer);
	VTK_CREATE(RenderWindow, main_window);
	main_window->SetRenderer(main_renderer);
	VTK_CREATE(RenderWindowInteractor, main_interactor);

	VTK_CREATE(Renderer, x_renderer);
	VTK_CREATE(RenderWindow, x_window);
	x_window->SetRenderer(x_renderer);
	VTK_CREATE(RenderWindowInteractor, x_interactor);

	VTK_CREATE(Renderer, y_renderer);
	VTK_CREATE(RenderWindow, y_window);
	y_window->SetRenderer(y_renderer);
	VTK_CREATE(RenderWindowInteractor, y_interactor);

	VTK_CREATE(Renderer, z_renderer);
	VTK_CREATE(RenderWindow, z_window);
	z_window->SetRenderer(z_renderer);
	VTK_CREATE(RenderWindowInteractor, z_interactor);

	this->ui->vtk_main->SetRenderWindow(main_window);
	this->ui->vtk_x->SetRenderWindow(x_window);
	this->ui->vtk_y->SetRenderWindow(y_window);
	this->ui->vtk_z->SetRenderWindow(z_window);

	this->ui->show0_checkbox->setCheckState(Qt::Unchecked);
	this->ui->show1_checkbox->setCheckState(Qt::Unchecked);

	this->ui->label_edit->setText(QString("N/A"));
	this->ui->bounds_text->setText(QString("N/A"));

	this->ui->value0_spinbox->setValue(0);
	this->ui->value1_spinbox->setValue(0);

	this->ui->min_spinbox->setValue(0);
	this->ui->max_spinbox->setValue(100.);

	this->ui->ID_spinbox->setValue(0);

	// connections
	connect(this->ui->value0_spinbox, SIGNAL(valueChanged(double)),
            this, SLOT(slot_set_value0(double)));
	connect(this->ui->value1_spinbox, SIGNAL(valueChanged(double)),
            this, SLOT(slot_set_value1(double)));
	connect(this->ui->min_spinbox, SIGNAL(valueChanged(double)),
            this, SLOT(slot_set_min(double)));
	connect(this->ui->max_spinbox, SIGNAL(valueChanged(double)),
            this, SLOT(slot_set_max(double)));
	connect(this->ui->ID_spinbox, SIGNAL(valueChanged(int)),
            this, SLOT(slot_set_ID(int)));
	connect(this->ui->compute_push, SIGNAL(clicked()),
            this, SLOT(slot_compute()));
	connect(this->ui->front_push, SIGNAL(stateChanged(int)),
            this, SLOT(slot_front_view(int)));
	connect(this->ui->back_push, SIGNAL(stateChanged(int)),
            this, SLOT(slot_back_view(int)));
	connect(this->ui->top_push, SIGNAL(stateChanged(int)),
            this, SLOT(slot_top_view(int)));
	connect(this->ui->bottom_push, SIGNAL(stateChanged(int)),
            this, SLOT(slot_bottom_view(int)));
	connect(this->ui->left_push, SIGNAL(stateChanged(int)),
            this, SLOT(slot_left_view(int)));
	connect(this->ui->right_push, SIGNAL(stateChanged(int)),
            this, SLOT(slot_right_view(int)));



}
