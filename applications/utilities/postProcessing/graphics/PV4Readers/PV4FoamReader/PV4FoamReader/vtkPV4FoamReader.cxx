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
#include "vtkPV4FoamReader.h"

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
#include "vtkPV4Foam.H"

#undef EXPERIMENTAL_TIME_CACHING

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

vtkStandardNewMacro(vtkPV4FoamReader);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vtkPV4FoamReader::vtkPV4FoamReader()
{
    Debug = 0;
    vtkDebugMacro(<<"Constructor");

    SetNumberOfInputPorts(0);

    FileName  = NULL;
    foamData_ = NULL;

    output0_  = NULL;

#ifdef VTKPV4FOAM_DUALPORT
    // Add second output for the Lagrangian
    this->SetNumberOfOutputPorts(2);
    vtkMultiBlockDataSet *lagrangian = vtkMultiBlockDataSet::New();
    lagrangian->ReleaseData();

    this->GetExecutive()->SetOutputData(1, lagrangian);
    lagrangian->Delete();
#endif

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

    UpdateGUI = 0;

    PartSelection = vtkDataArraySelection::New();
    VolFieldSelection = vtkDataArraySelection::New();
    PointFieldSelection = vtkDataArraySelection::New();
    LagrangianFieldSelection = vtkDataArraySelection::New();

    // Setup the selection callback to modify this object when an array
    // selection is changed.
    SelectionObserver = vtkCallbackCommand::New();
    SelectionObserver->SetCallback
    (
        &vtkPV4FoamReader::SelectionModifiedCallback
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

vtkPV4FoamReader::~vtkPV4FoamReader()
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
int vtkPV4FoamReader::RequestInformation
(
    vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector
)
{
    vtkDebugMacro(<<"RequestInformation");

    if (Foam::vtkPV4Foam::debug)
    {
        cout<<"REQUEST_INFORMATION\n";
    }

    if (!FileName)
    {
        vtkErrorMacro("FileName has to be specified!");
        return 0;
    }

    int nInfo = outputVector->GetNumberOfInformationObjects();

    if (Foam::vtkPV4Foam::debug)
    {
        cout<<"RequestInformation with " << nInfo << " item(s)\n";
        for (int infoI = 0; infoI < nInfo; ++infoI)
        {
            outputVector->GetInformationObject(infoI)->Print(cout);
        }
    }

    if (!foamData_)
    {
        foamData_ = new Foam::vtkPV4Foam(FileName, this);
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
        foamData_ = NULL;
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

        if (Foam::vtkPV4Foam::debug > 1)
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
int vtkPV4FoamReader::RequestData
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

    if (Foam::vtkPV4Foam::debug)
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
    // to the vtkPV4Foam::setTime() method
    for (int infoI = 0; infoI < nInfo; ++infoI)
    {
        vtkInformation *outInfo = outputVector->GetInformationObject(infoI);

        if
        (
            outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())
         && outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS()) > 0
        )
        {
            requestTime[nRequestTime++] =
                outInfo->Get
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

    if (Foam::vtkPV4Foam::debug)
    {
        cout<< "update output with "
            << output->GetNumberOfBlocks() << " blocks\n";
    }


#ifdef EXPERIMENTAL_TIME_CACHING
    bool needsUpdate = false;

    if (!output0_)
    {
        output0_ = vtkMultiBlockDataSet::New();
        needsUpdate = true;
    }

    // This experimental bit of code seems to work for the geometry,
    // but trashes the fields and still triggers the GeometryFilter
    if (needsUpdate)
    {
        foamData_->Update(output);
        output0_->ShallowCopy(output);
    }
    else
    {
        output->ShallowCopy(output0_);
    }

    if (Foam::vtkPV4Foam::debug)
    {
        if (needsUpdate)
        {
            cout<< "full UPDATE ---------\n";
        }
        else
        {
            cout<< "cached UPDATE ---------\n";
        }

        cout<< "UPDATED output: ";
        output->Print(cout);

        cout<< "UPDATED output0_: ";
        output0_->Print(cout);
    }

#else

#ifdef VTKPV4FOAM_DUALPORT
    foamData_->Update
    (
        output,
        vtkMultiBlockDataSet::SafeDownCast
        (
            outputVector->GetInformationObject(1)->Get
            (
                vtkMultiBlockDataSet::DATA_OBJECT()
            )
        );
    );
#else
    foamData_->Update(output, output);
#endif

    updatePatchNamesView(ShowPatchNames);

#endif

    // Do any cleanup on the OpenFOAM side
    foamData_->CleanUp();

    return 1;
}


void vtkPV4FoamReader::SetRefresh(int val)
{
    Modified();
}


void vtkPV4FoamReader::SetIncludeSets(int val)
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


void vtkPV4FoamReader::SetIncludeZones(int val)
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


void vtkPV4FoamReader::SetShowPatchNames(int val)
{
    if (ShowPatchNames != val)
    {
        ShowPatchNames = val;
        updatePatchNamesView(ShowPatchNames);
    }
}


void vtkPV4FoamReader::SetShowGroupsOnly(int val)
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


void vtkPV4FoamReader::updatePatchNamesView(const bool show)
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


void vtkPV4FoamReader::PrintSelf(ostream& os, vtkIndent indent)
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


int vtkPV4FoamReader::GetTimeStep()
{
    return foamData_ ? foamData_->timeIndex() : -1;
}


// ----------------------------------------------------------------------
// Parts selection list control

vtkDataArraySelection* vtkPV4FoamReader::GetPartSelection()
{
    vtkDebugMacro(<<"GetPartSelection");
    return PartSelection;
}


int vtkPV4FoamReader::GetNumberOfPartArrays()
{
    vtkDebugMacro(<<"GetNumberOfPartArrays");
    return PartSelection->GetNumberOfArrays();
}


const char* vtkPV4FoamReader::GetPartArrayName(int index)
{
    vtkDebugMacro(<<"GetPartArrayName");
    return PartSelection->GetArrayName(index);
}


int vtkPV4FoamReader::GetPartArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetPartArrayStatus");
    return PartSelection->ArrayIsEnabled(name);
}


void vtkPV4FoamReader::SetPartArrayStatus(const char* name, int status)
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

vtkDataArraySelection* vtkPV4FoamReader::GetVolFieldSelection()
{
    vtkDebugMacro(<<"GetVolFieldSelection");
    return VolFieldSelection;
}


int vtkPV4FoamReader::GetNumberOfVolFieldArrays()
{
    vtkDebugMacro(<<"GetNumberOfVolFieldArrays");
    return VolFieldSelection->GetNumberOfArrays();
}


const char* vtkPV4FoamReader::GetVolFieldArrayName(int index)
{
    vtkDebugMacro(<<"GetVolFieldArrayName");
    return VolFieldSelection->GetArrayName(index);
}


int vtkPV4FoamReader::GetVolFieldArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetVolFieldArrayStatus");
    return VolFieldSelection->ArrayIsEnabled(name);
}


void vtkPV4FoamReader::SetVolFieldArrayStatus(const char* name, int status)
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

vtkDataArraySelection* vtkPV4FoamReader::GetPointFieldSelection()
{
    vtkDebugMacro(<<"GetPointFieldSelection");
    return PointFieldSelection;
}


int vtkPV4FoamReader::GetNumberOfPointFieldArrays()
{
    vtkDebugMacro(<<"GetNumberOfPointFieldArrays");
    return PointFieldSelection->GetNumberOfArrays();
}


const char* vtkPV4FoamReader::GetPointFieldArrayName(int index)
{
    vtkDebugMacro(<<"GetPointFieldArrayName");
    return PointFieldSelection->GetArrayName(index);
}


int vtkPV4FoamReader::GetPointFieldArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetPointFieldArrayStatus");
    return PointFieldSelection->ArrayIsEnabled(name);
}


void vtkPV4FoamReader::SetPointFieldArrayStatus(const char* name, int status)
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

vtkDataArraySelection* vtkPV4FoamReader::GetLagrangianFieldSelection()
{
    vtkDebugMacro(<<"GetLagrangianFieldSelection");
    return LagrangianFieldSelection;
}


int vtkPV4FoamReader::GetNumberOfLagrangianFieldArrays()
{
    vtkDebugMacro(<<"GetNumberOfLagrangianFieldArrays");
    return LagrangianFieldSelection->GetNumberOfArrays();
}


const char* vtkPV4FoamReader::GetLagrangianFieldArrayName(int index)
{
    vtkDebugMacro(<<"GetLagrangianFieldArrayName");
    return LagrangianFieldSelection->GetArrayName(index);
}


int vtkPV4FoamReader::GetLagrangianFieldArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetLagrangianFieldArrayStatus");
    return LagrangianFieldSelection->ArrayIsEnabled(name);
}


void vtkPV4FoamReader::SetLagrangianFieldArrayStatus
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

void vtkPV4FoamReader::SelectionModifiedCallback
(
    vtkObject*,
    unsigned long,
    void* clientdata,
    void*
)
{
    static_cast<vtkPV4FoamReader*>(clientdata)->SelectionModified();
}


void vtkPV4FoamReader::SelectionModified()
{
    vtkDebugMacro(<<"SelectionModified");
    Modified();
}


int vtkPV4FoamReader::FillOutputPortInformation
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
