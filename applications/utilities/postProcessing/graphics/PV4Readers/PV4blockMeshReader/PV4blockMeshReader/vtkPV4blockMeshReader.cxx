/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
#include "vtkPV4blockMeshReader.h"

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
#include "vtkPV4blockMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

vtkStandardNewMacro(vtkPV4blockMeshReader);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vtkPV4blockMeshReader::vtkPV4blockMeshReader()
{
    Debug = 0;
    vtkDebugMacro(<<"Constructor");

    SetNumberOfInputPorts(0);

    FileName  = NULL;
    foamData_ = NULL;

    ShowPointNumbers = 1;
    UpdateGUI = 0;

    BlockSelection = vtkDataArraySelection::New();
    CurvedEdgesSelection = vtkDataArraySelection::New();

    // Setup the selection callback to modify this object when an array
    // selection is changed.
    SelectionObserver = vtkCallbackCommand::New();
    SelectionObserver->SetCallback
    (
        &vtkPV4blockMeshReader::SelectionModifiedCallback
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

vtkPV4blockMeshReader::~vtkPV4blockMeshReader()
{
    vtkDebugMacro(<<"Deconstructor");

    if (foamData_)
    {
        // remove point numbers
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

// Do everything except set the output info
int vtkPV4blockMeshReader::RequestInformation
(
    vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector
)
{
    vtkDebugMacro(<<"RequestInformation");

    if (Foam::vtkPV4blockMesh::debug)
    {
        cout<<"REQUEST_INFORMATION\n";
    }

    if (!FileName)
    {
        vtkErrorMacro("FileName has to be specified!");
        return 0;
    }

    int nInfo = outputVector->GetNumberOfInformationObjects();

    if (Foam::vtkPV4blockMesh::debug)
    {
        cout<<"RequestInformation with " << nInfo << " item(s)\n";
        for (int infoI = 0; infoI < nInfo; ++infoI)
        {
            outputVector->GetInformationObject(infoI)->Print(cout);
        }
    }

    if (!foamData_)
    {
        foamData_ = new Foam::vtkPV4blockMesh(FileName, this);
    }
    else
    {
        foamData_->updateInfo();
    }

    // might need some other type of error handling

//    {
//        vtkErrorMacro("could not find valid OpenFOAM blockMesh");
//
//        // delete foamData and flag it as fatal error
//        delete foamData_;
//        foamData_ = NULL;
//        return 0;
//    }


    return 1;
}


// Set the output info
int vtkPV4blockMeshReader::RequestData
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

    // catch previous error
    if (!foamData_)
    {
        vtkErrorMacro("Reader failed - perhaps no mesh?");
        return 0;
    }

    int nInfo = outputVector->GetNumberOfInformationObjects();

    if (Foam::vtkPV4blockMesh::debug)
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

    if (Foam::vtkPV4blockMesh::debug)
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



void vtkPV4blockMeshReader::SetShowPointNumbers(const int val)
{
    if (ShowPointNumbers != val)
    {
        ShowPointNumbers = val;
        updatePointNumbersView(ShowPointNumbers);
    }
}


void vtkPV4blockMeshReader::updatePointNumbersView(const bool show)
{
    pqApplicationCore* appCore = pqApplicationCore::instance();

    // need to check this, since our destructor calls this
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

    // use refresh here?
}


void vtkPV4blockMeshReader::PrintSelf(ostream& os, vtkIndent indent)
{
    vtkDebugMacro(<<"PrintSelf");

    this->Superclass::PrintSelf(os,indent);
    os  << indent << "File name: "
        << (this->FileName ? this->FileName : "(none)") << "\n";

    foamData_->PrintSelf(os, indent);
}


// ----------------------------------------------------------------------
// Block selection list control

vtkDataArraySelection* vtkPV4blockMeshReader::GetBlockSelection()
{
    vtkDebugMacro(<<"GetBlockSelection");
    return BlockSelection;
}


int vtkPV4blockMeshReader::GetNumberOfBlockArrays()
{
    vtkDebugMacro(<<"GetNumberOfBlockArrays");
    return BlockSelection->GetNumberOfArrays();
}


const char* vtkPV4blockMeshReader::GetBlockArrayName(int index)
{
    vtkDebugMacro(<<"GetBlockArrayName");
    return BlockSelection->GetArrayName(index);
}


int vtkPV4blockMeshReader::GetBlockArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetBlockArrayStatus");
    return BlockSelection->ArrayIsEnabled(name);
}


void vtkPV4blockMeshReader::SetBlockArrayStatus
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

vtkDataArraySelection* vtkPV4blockMeshReader::GetCurvedEdgesSelection()
{
    vtkDebugMacro(<<"GetCurvedEdgesSelection");
    return CurvedEdgesSelection;
}


int vtkPV4blockMeshReader::GetNumberOfCurvedEdgesArrays()
{
    vtkDebugMacro(<<"GetNumberOfCurvedEdgesArrays");
    return CurvedEdgesSelection->GetNumberOfArrays();
}


const char* vtkPV4blockMeshReader::GetCurvedEdgesArrayName(int index)
{
    vtkDebugMacro(<<"GetCurvedEdgesArrayName");
    return CurvedEdgesSelection->GetArrayName(index);
}


int vtkPV4blockMeshReader::GetCurvedEdgesArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetCurvedEdgesArrayStatus");
    return CurvedEdgesSelection->ArrayIsEnabled(name);
}


void vtkPV4blockMeshReader::SetCurvedEdgesArrayStatus
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

void vtkPV4blockMeshReader::SelectionModifiedCallback
(
    vtkObject*,
    unsigned long,
    void* clientdata,
    void*
)
{
    static_cast<vtkPV4blockMeshReader*>(clientdata)->Modified();
}


int vtkPV4blockMeshReader::FillOutputPortInformation
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
