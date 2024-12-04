#ifndef _CORN_MRI_MAIN_WINDOW_HPP_
#define _CORN_MRI_MAIN_WINDOW_HPP_

#include <QMainWindow>
#include <ui_main_window.h>
#include <vtkSmartPointer.h>

class vtkActor;
class vtkRenderWindow;
class vtkObject;
class vtkDataSetMapper;

class QAction;

// Forward Qt class declarations
class Ui_main_window;

class corn_mri_main_window : public QWidget
{    
    Q_OBJECT
public:
    corn_mri_main_window();
    ~corn_mri_main_window();

private slots:
	void slot_set_value0(double);
	void slot_set_value1(double);
	void slot_set_min(double);
	void slot_set_max(double);
	void slot_show_value0(double);
	void slot_show_value1(double);
	void slot_set_label(Qstring);
	void slot_set_id(QString);
	void slot_front_view();
	void slot_back_view();
	void slot_top_push();
	void slot_bottom_view();
	void slot_left_view();
	void slot_right_view();

signals:

private:
	
public:
	Ui_main_window *ui;
};

#endif // _CORN_MRI_MAIN_WINDOW_HPP_