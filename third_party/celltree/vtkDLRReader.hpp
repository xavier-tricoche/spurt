#ifndef __vtkDLRReader_h
#define __vtkDLRReader_h

#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkStdString.h"

class vtkDataArraySelection;
class vtkStdString;

class VTK_IO_EXPORT vtkDLRReader : public vtkUnstructuredGridAlgorithm
{
public:
  
    static vtkDLRReader *New();
    vtkTypeRevisionMacro(vtkDLRReader, vtkUnstructuredGridAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    // Description:
    // Specify the name of the DLR file to read
    vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

    // Description:
    // Get the reader's output
    vtkUnstructuredGrid *GetOutput();
    vtkUnstructuredGrid *GetOutput(int index);

    // Description:
    // The following methods allow selective reading of variable.
    // By default, ALL data fields on the nodes are read, but this can
    // be modified.
    int GetNumberOfVariables();
    const char* GetVariableName(int index);
    int GetVariableStatus(const char* name);
    void SetVariableStatus(const char* name, int status);
    void DisableAllVariables();
    void EnableAllVariables();

protected:

    vtkDLRReader();
    ~vtkDLRReader();

    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    char* FileName;


    vtkStdString GridFileName;
    vtkStdString PvalFileName;

    vtkDataArraySelection* VariableSelection;

    vtkStdString *VariableName;
    vtkIdType *ComponentNumber; 
    vtkIdType  NumberOfVariables;

private:

    vtkDLRReader(const vtkDLRReader&);  // Not implemented.
    void operator=(const vtkDLRReader&);  // Not implemented.

    void ReadFile(vtkUnstructuredGrid *output);

    bool ReadCells( int, vtkUnstructuredGrid* );
    bool ReadPoints( int, vtkUnstructuredGrid* );
    bool ReadVelocity( int, vtkUnstructuredGrid* );
};

#endif

