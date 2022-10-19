#ifndef _vtkCellTree_h
#define _vtkCellTree_h

#include <vtkAbstractCellLocator.h>

class celltree;

class vtkCellTree : public vtkAbstractCellLocator 
{
public:

    vtkTypeRevisionMacro(vtkCellTree,vtkAbstractCellLocator);
    void PrintSelf(ostream& os, vtkIndent indent);
    static vtkCellTree *New();

    using vtkAbstractCellLocator::IntersectWithLine;
    using vtkAbstractCellLocator::FindClosestPoint;
    using vtkAbstractCellLocator::FindClosestPointWithinRadius;

    void FreeSearchStructure();
    void BuildLocator();

    void GenerateRepresentation( int level, vtkPolyData *pd );
    void GenerateRepresentationLeafs( vtkPolyData *pd );

    virtual int IntersectWithLine( double p1[3], double p2[3], 
                                   double tol, double& t, double x[3],
                                   double pcoords[3], int &subId )
    { 
        return this->Superclass::IntersectWithLine( p1, p2, tol, t, x, pcoords, subId ); 
    }

    virtual int IntersectWithLine( double p1[3], double p2[3], 
                                   double tol, double &t, double x[3],
                                   double pcoords[3], int &subId, vtkIdType &cellId );

    virtual int IntersectWithLine( double p1[3], double p2[3], 
                                   double tol, double &t, double x[3], 
                                   double pcoords[3], int &subId, vtkIdType &cellId, 
                                   vtkGenericCell *cell );

    // Description:
    // Take the passed line segment and intersect it with the data set.
    // This method assumes that the data set is a vtkPolyData that describes
    // a closed surface, and the intersection points that are returned in
    // 'points' alternate between entrance points and exit points.
    // The return value of the function is 0 if no intersections were found,
    // -1 if point 'a0' lies inside the closed surface, or +1 if point 'a0'
    // lies outside the closed surface.
    // Either 'points' or 'cellIds' can be set to NULL if you don't want
    // to receive that information. This method is currently only implemented
    // in vtkOBBTree
    virtual int IntersectWithLine( const double p1[3], const double p2[3],
                                   vtkPoints *points, vtkIdList *cellIds )
    { 
        return this->Superclass::IntersectWithLine(p1, p2, points, cellIds); 
    }

    // Description:
    // Returns the Id of the cell containing the point, 
    // returns -1 if no cell found. This interface uses a tolerance of zero
    virtual vtkIdType FindCell(double x[3])
    { 
        return this->Superclass::FindCell(x); 
    }

    // Description:
    // Test a point to find if it is inside a cell. Returns the cellId if inside
    // or -1 if not.
    virtual vtkIdType FindCell( double x[3], double tol2, vtkGenericCell *GenCell, 
                                double pcoords[3], double *weights );

    virtual unsigned long GetMemorySize();

protected:

    vtkCellTree();
    ~vtkCellTree();

    celltree* ct;

    void BuildLocatorIfNeeded();
    void ForceBuildLocator();
    void BuildLocatorInternal();

private:

    vtkCellTree(const vtkCellTree&);
    void operator=(const vtkCellTree&);
};

#endif


