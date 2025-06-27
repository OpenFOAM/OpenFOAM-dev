/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    vtkPVFoamReader

Description
    reads a dataset in OpenFOAM format

    vtkPVblockMeshReader creates an multiblock dataset.
    It uses the OpenFOAM infrastructure (fvMesh, etc) to handle mesh and
    field data.

SourceFiles
    vtkPVblockMeshReader.cxx

\*---------------------------------------------------------------------------*/
#ifndef vtkPVFoamReader_h
#define vtkPVFoamReader_h

// VTK includes
#include "vtkMultiBlockDataSetAlgorithm.h"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

// VTK forward declarations
class vtkDataArraySelection;
class vtkCallbackCommand;

// OpenFOAM forward declarations
namespace Foam
{
    class vtkPVFoam;
}


/*---------------------------------------------------------------------------*\
                       Class vtkPVFoamReader Declaration
\*---------------------------------------------------------------------------*/

class vtkPVFoamReader
:
    public vtkMultiBlockDataSetAlgorithm
{
public:
    vtkTypeMacro(vtkPVFoamReader, vtkMultiBlockDataSetAlgorithm);

    void PrintSelf(ostream&, vtkIndent);

    static vtkPVFoamReader* New();

    // Description:
    // Set/Get the filename.
    vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

    // Description:
    // OpenFOAM read the decomposed case from the processor directories
    vtkSetMacro(DecomposedCase, int);
    vtkGetMacro(DecomposedCase, int);

    // Description:
    // OpenFOAM mesh caching control
    vtkSetMacro(CacheMesh, int);
    vtkGetMacro(CacheMesh, int);

    // Description:
    // OpenFOAM refresh times/fields
    virtual void SetRefresh();

    // Description:
    // OpenFOAM skip/include the 0/ time directory
    vtkSetMacro(SkipZeroTime, int);
    vtkGetMacro(SkipZeroTime, int);

    // Description:
    // OpenFOAM extrapolate internal values onto the patches
    vtkSetMacro(ExtrapolatePatches, int);
    vtkGetMacro(ExtrapolatePatches, int);

    // Description:
    // OpenFOAM use vtkPolyhedron instead of decomposing polyhedra
    vtkSetMacro(UseVTKPolyhedron, int);
    vtkGetMacro(UseVTKPolyhedron, int);

    // Description:
    // OpenFOAM read sets control
    virtual void SetIncludeSets(int);
    vtkGetMacro(IncludeSets, int);

    // Description:
    // OpenFOAM read zones control
    virtual void SetIncludeZones(int);
    vtkGetMacro(IncludeZones, int);

    // Description:
    // OpenFOAM display patch names control
    virtual void SetShowPatchNames(int);
    vtkGetMacro(ShowPatchNames, int);

    // Description:
    // OpenFOAM display patchGroups
    virtual void SetShowGroupsOnly(int);
    vtkGetMacro(ShowGroupsOnly, int);

    // Description:
    // OpenFOAM volField interpolation
    vtkSetMacro(InterpolateVolFields, int);
    vtkGetMacro(InterpolateVolFields, int);

    // Description:
    // Get the current timestep
    int GetTimeStep();

    // Description:
    // Parts selection list control
    virtual vtkDataArraySelection* GetPartSelection();
    int GetNumberOfPartArrays();
    int GetPartArrayStatus(const char* name);
    void SetPartArrayStatus(const char* name, int status);
    const char* GetPartArrayName(int index);

    // Description:
    // Field selection list control
    virtual vtkDataArraySelection* GetFieldSelection();
    int GetNumberOfFieldArrays();
    int GetFieldArrayStatus(const char* name);
    void SetFieldArrayStatus(const char* name, int status);
    const char* GetFieldArrayName(int index);

    // Description:
    // lagrangianField selection list control
    virtual vtkDataArraySelection* GetlagrangianFieldSelection();
    int GetNumberOflagrangianFieldArrays();
    int GetlagrangianFieldArrayStatus(const char* name);
    void SetlagrangianFieldArrayStatus(const char* name, int status);
    const char* GetlagrangianFieldArrayName(int index);

    // Description:
    // LagrangianField selection list control
    virtual vtkDataArraySelection* GetLagrangianFieldSelection();
    int GetNumberOfLagrangianFieldArrays();
    int GetLagrangianFieldArrayStatus(const char* name);
    void SetLagrangianFieldArrayStatus(const char* name, int status);
    const char* GetLagrangianFieldArrayName(int index);

    // Description:
    // Callback registered with the SelectionModifiedObserver
    // for all the selection lists
    static void SelectionModifiedCallback
    (
        vtkObject* caller,
        unsigned long eid,
        void* clientdata,
        void* calldata
    );

    // Description:
    // Callback registered with the AnimationPropertyModifiedObserver
    // for changes to animation properties
    static void AnimationPropertyModifiedCallback
    (
        vtkObject* caller,
        unsigned long eid,
        void* clientdata,
        void* calldata
    );


protected:

    //- Construct null
    vtkPVFoamReader();

    //- Disallow default bitwise copy construct
    vtkPVFoamReader(const vtkPVFoamReader&) = delete;

    //- Destructor
    ~vtkPVFoamReader();

    //- Return information about mesh, times, etc without loading anything
    virtual int RequestInformation
    (
        vtkInformation*,
        vtkInformationVector**,
        vtkInformationVector*
    );

    //- Get the mesh/fields for a particular time
    virtual int RequestData
    (
        vtkInformation*,
        vtkInformationVector**,
        vtkInformationVector*
    );

    //- Fill in additional port information
    virtual int FillOutputPortInformation(int, vtkInformation*);

    //- The observer to modify this object when array selections are modified
    vtkCallbackCommand* SelectionModifiedObserver;

    //- The observer to set the play-mode when animation properties are modified
    vtkCallbackCommand* AnimationPropertyModifiedObserver;

    //- The file name for this case
    char* FileName;

    //- Disallow default bitwise assignment
    void operator=(const vtkPVFoamReader&) = delete;


private:

    int Refresh;
    int DecomposedCase;
    int CacheMesh;
    int SkipZeroTime;

    int ExtrapolatePatches;
    int UseVTKPolyhedron;
    int IncludeSets;
    int IncludeZones;
    int ShowPatchNames;
    int ShowGroupsOnly;
    int InterpolateVolFields;

    bool First;

    vtkDataArraySelection* PartSelection;
    vtkDataArraySelection* FieldSelection;
    vtkDataArraySelection* lagrangianFieldSelection;
    vtkDataArraySelection* LagrangianFieldSelection;

    // BTX
    Foam::vtkPVFoam* foamData_;
    // ETX

    //- Return information about the available time steps
    int requestTimeSteps(vtkInformationVector*);

    //- Add/remove patch names to/from the view
    void updatePatchNamesView(const bool show);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
