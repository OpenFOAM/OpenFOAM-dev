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

\*---------------------------------------------------------------------------*/
#include "vtkPVblockMeshReader.h"

#include "pqApplicationCore.h"
#include "pqRenderView.h"
#include "pqServerManagerModel.h"

// VTK includes
#include "vtkCallbackCommand.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkObjectFactory.h"
#include "vtkSMRenderViewProxy.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"

// OpenFOAM includes
#include "vtkPVblockMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

vtkStandardNewMacro(vtkPVblockMeshReader);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vtkPVblockMeshReader::vtkPVblockMeshReader()
{
    Debug = 0;
    vtkDebugMacro(<<"Constructor");

    SetNumberOfInputPorts(0);

    FileName  = nullptr;
    foamData_ = nullptr;

    ShowPointNumbers = 1;

    BlockSelection = vtkDataArraySelection::New();
    CurvedEdgesSelection = vtkDataArraySelection::New();

    // Setup the selection callback to modify this object when an array
    // selection is changed.
    SelectionObserver = vtkCallbackCommand::New();
    SelectionObserver->SetCallback
    (
        &vtkPVblockMeshReader::SelectionModifiedCallback
    );
    SelectionObserver->SetClientData(this);


    BlockSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );

    CurvedEdgesSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

vtkPVblockMeshReader::~vtkPVblockMeshReader()
{
    vtkDebugMacro(<<"Deconstructor");

    if (foamData_)
    {
        // Remove point numbers
        updatePointNumbersView(false);
        delete foamData_;
    }

    if (FileName)
    {
        delete [] FileName;
    }

    BlockSelection->RemoveObserver(this->SelectionObserver);
    CurvedEdgesSelection->RemoveObserver(this->SelectionObserver);

    SelectionObserver->Delete();
    BlockSelection->Delete();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

int vtkPVblockMeshReader::RequestInformation
(
    vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector
)
{
    vtkDebugMacro(<<"RequestInformation");

    if (Foam::vtkPVblockMesh::debug)
    {
        cout<<"REQUEST_INFORMATION\n";
    }

    if (!FileName)
    {
        vtkErrorMacro("FileName has to be specified!");
        return 0;
    }

    int nInfo = outputVector->GetNumberOfInformationObjects();

    if (Foam::vtkPVblockMesh::debug)
    {
        cout<<"RequestInformation with " << nInfo << " item(s)\n";
        for (int infoI = 0; infoI < nInfo; ++infoI)
        {
            outputVector->GetInformationObject(infoI)->Print(cout);
        }
    }

    if (!foamData_)
    {
        foamData_ = new Foam::vtkPVblockMesh(FileName, this);
    }
    else
    {
        foamData_->updateInfo();
    }

    return 1;
}


int vtkPVblockMeshReader::RequestData
(
    vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector
)
{
    vtkDebugMacro(<<"RequestData");

    if (!FileName)
    {
        vtkErrorMacro("FileName has to be specified!");
        return 0;
    }

    // Catch previous error
    if (!foamData_)
    {
        vtkErrorMacro("Reader failed - perhaps no mesh?");
        return 0;
    }

    int nInfo = outputVector->GetNumberOfInformationObjects();

    if (Foam::vtkPVblockMesh::debug)
    {
        cout<<"RequestData with " << nInfo << " item(s)\n";
        for (int infoI = 0; infoI < nInfo; ++infoI)
        {
            outputVector->GetInformationObject(infoI)->Print(cout);
        }
    }

    vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast
    (
        outputVector->GetInformationObject(0)->Get
        (
            vtkMultiBlockDataSet::DATA_OBJECT()
        )
    );

    if (Foam::vtkPVblockMesh::debug)
    {
        cout<< "update output with "
            << output->GetNumberOfBlocks() << " blocks\n";
    }


    foamData_->Update(output);
    updatePointNumbersView(ShowPointNumbers);

    // Do any cleanup on the OpenFOAM side
    foamData_->CleanUp();

    return 1;
}


void vtkPVblockMeshReader::SetRefresh()
{
    // Delete the current blockMesh to force re-read and update
    if (foamData_)
    {
        updatePointNumbersView(false);
        delete foamData_;
        foamData_ = 0;
    }

    Modified();
}


void vtkPVblockMeshReader::SetShowPointNumbers(const int val)
{
    if (ShowPointNumbers != val)
    {
        ShowPointNumbers = val;
        updatePointNumbersView(ShowPointNumbers);
    }
}


void vtkPVblockMeshReader::updatePointNumbersView(const bool show)
{
    pqApplicationCore* appCore = pqApplicationCore::instance();

    // Need to check this, since our destructor calls this
    if (!appCore)
    {
        return;
    }

    // Server manager model for querying items in the server manager
    pqServerManagerModel* smModel = appCore->getServerManagerModel();
    if (!smModel || !foamData_)
    {
        return;
    }


    // Get all the pqRenderView instances
    QList<pqRenderView*> renderViews = smModel->findItems<pqRenderView*>();
    for (int viewI=0; viewI<renderViews.size(); ++viewI)
    {
        foamData_->renderPointNumbers
        (
            renderViews[viewI]->getRenderViewProxy()->GetRenderer(),
            show
        );
    }

    // Use refresh here?
}


void vtkPVblockMeshReader::PrintSelf(ostream& os, vtkIndent indent)
{
    vtkDebugMacro(<<"PrintSelf");

    this->Superclass::PrintSelf(os,indent);
    os  << indent << "File name: "
        << (this->FileName ? this->FileName : "(none)") << "\n";

    foamData_->PrintSelf(os, indent);
}


// ----------------------------------------------------------------------
// Block selection list control

vtkDataArraySelection* vtkPVblockMeshReader::GetBlockSelection()
{
    vtkDebugMacro(<<"GetBlockSelection");
    return BlockSelection;
}


int vtkPVblockMeshReader::GetNumberOfBlockArrays()
{
    vtkDebugMacro(<<"GetNumberOfBlockArrays");
    return BlockSelection->GetNumberOfArrays();
}


const char* vtkPVblockMeshReader::GetBlockArrayName(int index)
{
    vtkDebugMacro(<<"GetBlockArrayName");
    return BlockSelection->GetArrayName(index);
}


int vtkPVblockMeshReader::GetBlockArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetBlockArrayStatus");
    return BlockSelection->ArrayIsEnabled(name);
}


void vtkPVblockMeshReader::SetBlockArrayStatus
(
    const char* name,
    int status
)
{
    vtkDebugMacro(<<"SetBlockArrayStatus");
    if (status)
    {
        BlockSelection->EnableArray(name);
    }
    else
    {
        BlockSelection->DisableArray(name);
    }
}


// ----------------------------------------------------------------------
// CurvedEdges selection list control

vtkDataArraySelection* vtkPVblockMeshReader::GetCurvedEdgesSelection()
{
    vtkDebugMacro(<<"GetCurvedEdgesSelection");
    return CurvedEdgesSelection;
}


int vtkPVblockMeshReader::GetNumberOfCurvedEdgesArrays()
{
    vtkDebugMacro(<<"GetNumberOfCurvedEdgesArrays");
    return CurvedEdgesSelection->GetNumberOfArrays();
}


const char* vtkPVblockMeshReader::GetCurvedEdgesArrayName(int index)
{
    vtkDebugMacro(<<"GetCurvedEdgesArrayName");
    return CurvedEdgesSelection->GetArrayName(index);
}


int vtkPVblockMeshReader::GetCurvedEdgesArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetCurvedEdgesArrayStatus");
    return CurvedEdgesSelection->ArrayIsEnabled(name);
}


void vtkPVblockMeshReader::SetCurvedEdgesArrayStatus
(
    const char* name,
    int status
)
{
    vtkDebugMacro(<<"SetCurvedEdgesArrayStatus");
    if (status)
    {
        CurvedEdgesSelection->EnableArray(name);
    }
    else
    {
        CurvedEdgesSelection->DisableArray(name);
    }
}


// ----------------------------------------------------------------------

void vtkPVblockMeshReader::SelectionModifiedCallback
(
    vtkObject*,
    unsigned long,
    void* clientdata,
    void*
)
{
    static_cast<vtkPVblockMeshReader*>(clientdata)->Modified();
}


int vtkPVblockMeshReader::FillOutputPortInformation
(
    int port,
    vtkInformation* info
)
{
    if (port == 0)
    {
        return this->Superclass::FillOutputPortInformation(port, info);
    }
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
}


// ************************************************************************* //
