/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    vtkPVblockMeshReader

Description
    reads a dataset in OpenFOAM bockMesh format

    vtkPVblockMeshReader creates an multiblock dataset.
    It uses the OpenFOAM infrastructure (blockMesh).

SourceFiles
    vtkPVblockMeshReader.cxx

\*---------------------------------------------------------------------------*/

#ifndef vtkPVblockMeshReader_h
#define vtkPVblockMeshReader_h

// VTK includes
#include "vtkMultiBlockDataSetAlgorithm.h"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

// VTK forward declarations
class vtkDataArraySelection;
class vtkCallbackCommand;

namespace Foam
{
    class vtkPVblockMesh;
}

/*---------------------------------------------------------------------------*\
                   Class vtkPVblockMeshReader Declaration
\*---------------------------------------------------------------------------*/

class vtkPVblockMeshReader
:
    public vtkMultiBlockDataSetAlgorithm
{
public:
    vtkTypeMacro(vtkPVblockMeshReader, vtkMultiBlockDataSetAlgorithm);
    void PrintSelf(ostream&, vtkIndent);

    static vtkPVblockMeshReader* New();

    // Description:
    // Set/Get the filename.
    vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

    // Description:
    // Display corner point labels
    virtual void SetShowPointNumbers(int);
    vtkGetMacro(ShowPointNumbers, int);

    // Description:
    // Refresh blockMesh from changes to blockMeshDict
    virtual void SetRefresh();

    // Description:
    // Blocks selection list control
    vtkDataArraySelection* GetBlockSelection();
    int  GetNumberOfBlockArrays();
    int  GetBlockArrayStatus(const char*);
    void SetBlockArrayStatus(const char*, int status);
    const char* GetBlockArrayName(int index);

    // Description:
    // CurvedEdges selection list control
    vtkDataArraySelection* GetCurvedEdgesSelection();
    int  GetNumberOfCurvedEdgesArrays();
    int  GetCurvedEdgesArrayStatus(const char*);
    void SetCurvedEdgesArrayStatus(const char*, int status);
    const char* GetCurvedEdgesArrayName(int index);

    // Description:
    // Callback registered with the SelectionObserver
    // for all the selection lists
    static void SelectionModifiedCallback
    (
        vtkObject* caller,
        unsigned long eid,
        void* clientdata,
        void* calldata
    );


protected:

    //- Construct null
    vtkPVblockMeshReader();

    //- Destructor
    ~vtkPVblockMeshReader();

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

    // The observer to modify this object when array selections are modified
    vtkCallbackCommand* SelectionObserver;

    char* FileName;


private:

    //- Disallow default bitwise copy construct
    vtkPVblockMeshReader(const vtkPVblockMeshReader&);

    //- Disallow default bitwise assignment
    void operator=(const vtkPVblockMeshReader&);

    //- Add/remove point numbers to/from the view
    void updatePointNumbersView(const bool show);


    //- Show Point Numbers
    int ShowPointNumbers;

    vtkDataArraySelection* BlockSelection;

    vtkDataArraySelection* CurvedEdgesSelection;

    // BTX
    Foam::vtkPVblockMesh* foamData_;
    // ETX
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
