#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLightKit.h"
#include "vtkPolyDataNormals.h"
#include "vtkCommand.h"
#include "vtkInteractorObserver.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkTransform.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGReader.h"
#include "vtkPNGWriter.h"
#include "vtkGlyph3D.h"
#include "vtkSource.h"
#include "vtkSphereSource.h"
#include "vtkPoints.h"
#include "vtkFloatArray.h"
#include "vtkImageActor.h"
#include "vtkSmartPointer.h"
#include "vtkImageMapper3D.h"


#include <string>
#include <util/timer.hpp>
#include <sstream>
#include <math/fixed_vector.hpp>
#include <reconstruction/XY2LatLong.hpp>


void export_camera_info(vtkRenderer* ren)
{
    nvis::vec3 foc, pos, up;
    nvis::vec2 clip;
    ren->GetActiveCamera()->GetPosition(pos.begin());
    ren->GetActiveCamera()->GetFocalPoint(foc.begin());
    ren->GetActiveCamera()->GetViewUp(up.begin());
    ren->GetActiveCamera()->GetClippingRange(clip.begin());
    std::cout
            << "# Camera info:\n"
            << "position " << pos << '\n'
            << "focal_point " << foc << '\n'
            << "up_vector " << up << '\n'
            << "clipping_range " << clip << '\n';
}

// Callback for the interaction
class cameraCB : public vtkCommand {
public:
    static cameraCB* New() {
        return new cameraCB;
    }
    virtual void Execute(vtkObject* caller, unsigned long, void*) {
        if (true) {
            vtkRenderer* ren = reinterpret_cast<vtkRenderer*>(caller);
            export_camera_info(ren);
        }
    }
};

class timerCB : public vtkCommand {
public:
    static timerCB* New() {
        timerCB* cb = new timerCB;
        cb->TimerCount = 0;
        return cb;
    }
    
    virtual void Execute(vtkObject* caller, unsigned long eventId,
                         void* vtkNotUsed(callData)) {
        if (vtkCommand::TimerEvent == eventId) {
            ++this->TimerCount;
        }
        std::cout << this->TimerCount << std::endl;
        actor->SetPosition(this->TimerCount, this->TimerCount,0);
        
        vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::SafeDownCast(caller);
        iren->GetRenderWindow()->Render();
    }
    
private:
    int TimerCount;
public:
    vtkActor* actor;
};

void readCSV(std::fstream& f, int& year, std::string& degree, float& lo, float& la)
{
    std::string c;
    f >> la >> c >> lo >> c;
    for (int i=0 ; i<2 ; ++i) {
        c = ' ';
        while (c != ',') {
            f >> c;
        }
    }
    f >> year >> c >> degree;
}

void printUsageAndExit( const std::string& argv0, bool doExit = true )
{
    std::cerr
            << "Usage  : " << argv0 << " [options]\n"
            << "Options:\n"
            << "    -h   | --help                      Print this information\n"
            << "    -m   | --map <string>              Map file name (PNG)\n"
            << "    -z   | --zoom <float>              Zoom factor for camera\n"
            << "    -c   | --cities <string>           Cities file name (CSV)\n"
            << "    -i   | --input <string>            Date file name (CSV)\n";
            << std::endl;
            
    if (doExit) {
        exit(1);
    }
}

int main(int argc, char* argv[])
{

    std::string filename = argv[1];
    
    vtkSmartPointer<vtkPNGReader> pngin = vtkSmartPointer<vtkPNGReader>::New();
    pngin->SetDataSpacing(1,1,1);
    pngin->SetFileName(filename.c_str());
    vtkSmartPointer<vtkImageActor> imactor = vtkSmartPointer<vtkImageActor>::New();
    imactor->GetMapper()->SetInputConnection(pngin->GetOutputPort());
    
    vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<cameraCB> cb = vtkSmartPointer<cameraCB>::New();
    ren->AddObserver(vtkCommand::StartEvent, cb);
    
    vtkSmartPointer<vtkRenderWindow> renwin = vtkSmartPointer<vtkRenderWindow>::New();
    renwin->AddRenderer(ren);
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renwin);
    iren->Initialize();
    ren->AddActor(imactor);
    ren->SetBackground(0,0,0);
    renwin->SetSize(1260, 609);
    
    ren->ResetCamera();
    vtkCamera* camera = ren->GetActiveCamera();
    if (argc == 3) {
        camera->Zoom(atof(argv[2]));
    }
    renwin->Render();
    
    iren->Start();
}
