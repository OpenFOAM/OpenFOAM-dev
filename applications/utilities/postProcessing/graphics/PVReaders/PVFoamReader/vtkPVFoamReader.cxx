/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
#include "vtkPVFoamReader.h"

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
#include "vtkPVFoam.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

vtkStandardNewMacro(vtkPVFoamReader);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vtkPVFoamReader::vtkPVFoamReader()
{
    Debug = 0;
    vtkDebugMacro(<<"Constructor");

    SetNumberOfInputPorts(0);

    FileName  = nullptr;
    foamData_ = nullptr;

    output0_  = nullptr;

    TimeStepRange[0] = 0;
    TimeStepRange[1] = 0;

    CacheMesh = 1;
    Refresh = 0;

    SkipZeroTime = 0;
    ExtrapolatePatches = 0;
    UseVTKPolyhedron = 0;
    IncludeSets = 0;
    IncludeZones = 0;
    ShowPatchNames = 0;
    ShowGroupsOnly = 0;
    InterpolateVolFields = 1;

    PartSelection = vtkDataArraySelection::New();
    VolFieldSelection = vtkDataArraySelection::New();
    PointFieldSelection = vtkDataArraySelection::New();
    LagrangianFieldSelection = vtkDataArraySelection::New();

    // Setup the selection callback to modify this object when an array
    // selection is changed.
    SelectionObserver = vtkCallbackCommand::New();
    SelectionObserver->SetCallback
    (
        &vtkPVFoamReader::SelectionModifiedCallback
    );
    SelectionObserver->SetClientData(this);

    PartSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
    VolFieldSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
    PointFieldSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
    LagrangianFieldSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

vtkPVFoamReader::~vtkPVFoamReader()
{
    vtkDebugMacro(<<"Deconstructor");

    if (foamData_)
    {
        // remove patch names
        updatePatchNamesView(false);
        delete foamData_;
    }

    if (FileName)
    {
        delete [] FileName;
    }

    if (output0_)
    {
        output0_->Delete();
    }


    PartSelection->RemoveObserver(this->SelectionObserver);
    VolFieldSelection->RemoveObserver(this->SelectionObserver);
    PointFieldSelection->RemoveObserver(this->SelectionObserver);
    LagrangianFieldSelection->RemoveObserver(this->SelectionObserver);

    SelectionObserver->Delete();

    PartSelection->Delete();
    VolFieldSelection->Delete();
    PointFieldSelection->Delete();
    LagrangianFieldSelection->Delete();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Do everything except set the output info
int vtkPVFoamReader::RequestInformation
(
    vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector
)
{
    vtkDebugMacro(<<"RequestInformation");

    if (Foam::vtkPVFoam::debug)
    {
        cout<<"REQUEST_INFORMATION\n";
    }

    if (!FileName)
    {
        vtkErrorMacro("FileName has to be specified!");
        return 0;
    }

    int nInfo = outputVector->GetNumberOfInformationObjects();

    if (Foam::vtkPVFoam::debug)
    {
        cout<<"RequestInformation with " << nInfo << " item(s)\n";
        for (int infoI = 0; infoI < nInfo; ++infoI)
        {
            outputVector->GetInformationObject(infoI)->Print(cout);
        }
    }

    if (!foamData_)
    {
        foamData_ = new Foam::vtkPVFoam(FileName, this);
    }
    else
    {
        foamData_->updateInfo();
    }

    int nTimeSteps = 0;
    double* timeSteps = foamData_->findTimes(nTimeSteps);

    if (!nTimeSteps)
    {
        vtkErrorMacro("could not find valid OpenFOAM mesh");

        // delete foamData and flag it as fatal error
        delete foamData_;
        foamData_ = nullptr;
        return 0;
    }

    // set identical time steps for all ports
    for (int infoI = 0; infoI < nInfo; ++infoI)
    {
        outputVector->GetInformationObject(infoI)->Set
        (
            vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
            timeSteps,
            nTimeSteps
        );
    }

    if (nTimeSteps)
    {
        double timeRange[2];
        timeRange[0] = timeSteps[0];
        timeRange[1] = timeSteps[nTimeSteps-1];

        if (Foam::vtkPVFoam::debug > 1)
        {
            cout<<"nTimeSteps " << nTimeSteps << "\n"
                <<"timeRange " << timeRange[0] << " to " << timeRange[1]
                << "\n";

            for (int timeI = 0; timeI < nTimeSteps; ++timeI)
            {
                cout<< "step[" << timeI << "] = " << timeSteps[timeI] << "\n";
            }
        }

        for (int infoI = 0; infoI < nInfo; ++infoI)
        {
            outputVector->GetInformationObject(infoI)->Set
            (
                vtkStreamingDemandDrivenPipeline::TIME_RANGE(),
                timeRange,
                2
            );
        }
    }

    delete timeSteps;

    return 1;
}


// Set the output info
int vtkPVFoamReader::RequestData
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

    if (!foamData_)
    {
        vtkErrorMacro("Reader failed - perhaps no mesh?");
        return 0;
    }

    int nInfo = outputVector->GetNumberOfInformationObjects();

    if (Foam::vtkPVFoam::debug)
    {
        cout<<"RequestData with " << nInfo << " item(s)\n";
        for (int infoI = 0; infoI < nInfo; ++infoI)
        {
            outputVector->GetInformationObject(infoI)->Print(cout);
        }
    }

    // Get the requested time step.
    // We only support requests for a single time step

    int nRequestTime = 0;
    double requestTime[nInfo];

    // taking port0 as the lead for other outputs would be nice, but fails when
    // a filter is added - we need to check everything
    // but since PREVIOUS_UPDATE_TIME_STEPS() is protected, relay the logic
    // to the vtkPVFoam::setTime() method
    for (int infoI = 0; infoI < nInfo; ++infoI)
    {
        vtkInformation *outInfo = outputVector->GetInformationObject(infoI);

        if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
        {
            int nTimes =
                outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

            requestTime[nRequestTime++] =
                nTimes == 1
              ? outInfo->Get
                (
                    vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
                    0
                )
              : outInfo->Get
              (
                  vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()
              );
        }
    }

    if (nRequestTime)
    {
        foamData_->setTime(nRequestTime, requestTime);
    }

    vtkMultiBlockDataSet* output = vtkMultiBlockDataSet::SafeDownCast
    (
        outputVector->GetInformationObject(0)->Get
        (
            vtkMultiBlockDataSet::DATA_OBJECT()
        )
    );

    if (Foam::vtkPVFoam::debug)
    {
        cout<< "update output with "
            << output->GetNumberOfBlocks() << " blocks\n";
    }

    foamData_->Update(output, output);

    updatePatchNamesView(ShowPatchNames);

    // Do any cleanup on the OpenFOAM side
    foamData_->CleanUp();

    return 1;
}


void vtkPVFoamReader::SetRefresh()
{
    Modified();

    pqApplicationCore::instance()->render();
}


void vtkPVFoamReader::SetIncludeSets(int val)
{
    if (IncludeSets != val)
    {
        IncludeSets = val;
        if (foamData_)
        {
            foamData_->updateInfo();
        }
    }
}


void vtkPVFoamReader::SetIncludeZones(int val)
{
    if (IncludeZones != val)
    {
        IncludeZones = val;
        if (foamData_)
        {
            foamData_->updateInfo();
        }
    }
}


void vtkPVFoamReader::SetShowPatchNames(int val)
{
    if (ShowPatchNames != val)
    {
        ShowPatchNames = val;
        updatePatchNamesView(ShowPatchNames);
    }
}


void vtkPVFoamReader::SetShowGroupsOnly(int val)
{
    if (ShowGroupsOnly != val)
    {
        ShowGroupsOnly = val;
        if (foamData_)
        {
            foamData_->updateInfo();
        }
    }
}


void vtkPVFoamReader::updatePatchNamesView(const bool show)
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

    for (int viewI=0; viewI < renderViews.size(); ++viewI)
    {
        foamData_->renderPatchNames
        (
            renderViews[viewI]->getRenderViewProxy()->GetRenderer(),
            show
        );
    }

    // use refresh here?
}


void vtkPVFoamReader::PrintSelf(ostream& os, vtkIndent indent)
{
    vtkDebugMacro(<<"PrintSelf");

    this->Superclass::PrintSelf(os,indent);
    os  << indent << "File name: "
        << (this->FileName ? this->FileName : "(none)") << "\n";

    foamData_->PrintSelf(os, indent);

    os  << indent << "Time step range: "
        << this->TimeStepRange[0] << " - " << this->TimeStepRange[1] << "\n"
        << indent << "Time step: " << this->GetTimeStep() << endl;
}


int vtkPVFoamReader::GetTimeStep()
{
    return foamData_ ? foamData_->timeIndex() : -1;
}


// ----------------------------------------------------------------------
// Parts selection list control

vtkDataArraySelection* vtkPVFoamReader::GetPartSelection()
{
    vtkDebugMacro(<<"GetPartSelection");
    return PartSelection;
}


int vtkPVFoamReader::GetNumberOfPartArrays()
{
    vtkDebugMacro(<<"GetNumberOfPartArrays");
    return PartSelection->GetNumberOfArrays();
}


const char* vtkPVFoamReader::GetPartArrayName(int index)
{
    vtkDebugMacro(<<"GetPartArrayName");
    return PartSelection->GetArrayName(index);
}


int vtkPVFoamReader::GetPartArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetPartArrayStatus");
    return PartSelection->ArrayIsEnabled(name);
}


void vtkPVFoamReader::SetPartArrayStatus(const char* name, int status)
{
    vtkDebugMacro("Set mesh part \"" << name << "\" status to: " << status);

    if (status)
    {
        PartSelection->EnableArray(name);
    }
    else
    {
        PartSelection->DisableArray(name);
    }
}


// ----------------------------------------------------------------------
// volField selection list control

vtkDataArraySelection* vtkPVFoamReader::GetVolFieldSelection()
{
    vtkDebugMacro(<<"GetVolFieldSelection");
    return VolFieldSelection;
}


int vtkPVFoamReader::GetNumberOfVolFieldArrays()
{
    vtkDebugMacro(<<"GetNumberOfVolFieldArrays");
    return VolFieldSelection->GetNumberOfArrays();
}


const char* vtkPVFoamReader::GetVolFieldArrayName(int index)
{
    vtkDebugMacro(<<"GetVolFieldArrayName");
    return VolFieldSelection->GetArrayName(index);
}


int vtkPVFoamReader::GetVolFieldArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetVolFieldArrayStatus");
    return VolFieldSelection->ArrayIsEnabled(name);
}


void vtkPVFoamReader::SetVolFieldArrayStatus(const char* name, int status)
{
    vtkDebugMacro(<<"SetVolFieldArrayStatus");
    if (status)
    {
        VolFieldSelection->EnableArray(name);
    }
    else
    {
        VolFieldSelection->DisableArray(name);
    }
}


// ----------------------------------------------------------------------
// pointField selection list control

vtkDataArraySelection* vtkPVFoamReader::GetPointFieldSelection()
{
    vtkDebugMacro(<<"GetPointFieldSelection");
    return PointFieldSelection;
}


int vtkPVFoamReader::GetNumberOfPointFieldArrays()
{
    vtkDebugMacro(<<"GetNumberOfPointFieldArrays");
    return PointFieldSelection->GetNumberOfArrays();
}


const char* vtkPVFoamReader::GetPointFieldArrayName(int index)
{
    vtkDebugMacro(<<"GetPointFieldArrayName");
    return PointFieldSelection->GetArrayName(index);
}


int vtkPVFoamReader::GetPointFieldArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetPointFieldArrayStatus");
    return PointFieldSelection->ArrayIsEnabled(name);
}


void vtkPVFoamReader::SetPointFieldArrayStatus(const char* name, int status)
{
    vtkDebugMacro(<<"SetPointFieldArrayStatus");
    if (status)
    {
        PointFieldSelection->EnableArray(name);
    }
    else
    {
        PointFieldSelection->DisableArray(name);
    }
}


// ----------------------------------------------------------------------
// lagrangianField selection list control

vtkDataArraySelection* vtkPVFoamReader::GetLagrangianFieldSelection()
{
    vtkDebugMacro(<<"GetLagrangianFieldSelection");
    return LagrangianFieldSelection;
}


int vtkPVFoamReader::GetNumberOfLagrangianFieldArrays()
{
    vtkDebugMacro(<<"GetNumberOfLagrangianFieldArrays");
    return LagrangianFieldSelection->GetNumberOfArrays();
}


const char* vtkPVFoamReader::GetLagrangianFieldArrayName(int index)
{
    vtkDebugMacro(<<"GetLagrangianFieldArrayName");
    return LagrangianFieldSelection->GetArrayName(index);
}


int vtkPVFoamReader::GetLagrangianFieldArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetLagrangianFieldArrayStatus");
    return LagrangianFieldSelection->ArrayIsEnabled(name);
}


void vtkPVFoamReader::SetLagrangianFieldArrayStatus
(
    const char* name,
    int status
)
{
    vtkDebugMacro(<<"SetLagrangianFieldArrayStatus");
    if (status)
    {
        LagrangianFieldSelection->EnableArray(name);
    }
    else
    {
        LagrangianFieldSelection->DisableArray(name);
    }
}


// ----------------------------------------------------------------------

void vtkPVFoamReader::SelectionModifiedCallback
(
    vtkObject*,
    unsigned long,
    void* clientdata,
    void*
)
{
    static_cast<vtkPVFoamReader*>(clientdata)->SelectionModified();
}


void vtkPVFoamReader::SelectionModified()
{
    vtkDebugMacro(<<"SelectionModified");
    Modified();
}


int vtkPVFoamReader::FillOutputPortInformation
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
