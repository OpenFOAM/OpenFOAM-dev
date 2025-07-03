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

\*---------------------------------------------------------------------------*/

#include "vtkPVFoamReader.h"

#include "pqApplicationCore.h"
#include "pqRenderView.h"
#include "pqServerManagerModel.h"
#include "pqPVApplicationCore.h"
#include "pqAnimationManager.h"
#include "pqAnimationScene.h"

// VTK includes
#include "vtkCallbackCommand.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkObjectFactory.h"
#include "vtkSMRenderViewProxy.h"
#include "vtkSMAnimationSceneProxy.h"
#include "vtkSMPropertyHelper.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include "vtkCompositeAnimationPlayer.h"

// OpenFOAM includes
#include "vtkPVFoam.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

int vtkPVFoamReader::requestTimeSteps(vtkInformationVector* outputVector)
{
    int nTimeSteps = 0;
    double* timeSteps = foamData_->findTimes(First, nTimeSteps);

    if (!nTimeSteps) return 0;

    int nInfo = outputVector->GetNumberOfInformationObjects();

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
            cout<< "\nnTimeSteps " << nTimeSteps
                << "\ntimeRange " << timeRange[0] << " to " << timeRange[1]
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


void vtkPVFoamReader::updatePatchNamesView(const bool show)
{
    if (!foamData_) return;

    pqApplicationCore* appCore = pqApplicationCore::instance();
    if (!appCore) return;

    // Server manager model for querying items in the server manager
    pqServerManagerModel* smModel = appCore->getServerManagerModel();
    if (!smModel) return;

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


vtkSMProxy* vtkPVFoamReader::getActiveAnimationSceneProxy()
{
    pqPVApplicationCore* appCore = pqPVApplicationCore::instance();
    if (!appCore) return nullptr;

    pqAnimationManager* animationManager = appCore->animationManager();
    if (!animationManager) return nullptr;

    pqAnimationScene* animationScene = animationManager->getActiveScene();
    if (!animationScene) return nullptr;

    return animationScene->getProxy();
}


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

    DecomposedCase = 0;
    CacheMesh = 1;
    Refresh = 0;

    SkipZeroTime = 0;
    ExtrapolatePatches = 0;
    UseVTKPolyhedron = 0;
    IncludeSets = 0;
    IncludeZones = 1;
    ShowPatchNames = 0;
    ShowGroupsOnly = 0;
    InterpolateVolFields = 1;

    First = true;

    PartSelection = vtkDataArraySelection::New();
    FieldSelection = vtkDataArraySelection::New();
    lagrangianFieldSelection = vtkDataArraySelection::New();
    LagrangianFieldSelection = vtkDataArraySelection::New();

    // Setup the selection-modified callback to modify this object when an
    // array selection is changed
    SelectionModifiedObserver = vtkCallbackCommand::New();
    SelectionModifiedObserver->SetCallback
    (
        &vtkPVFoamReader::SelectionModifiedCallback
    );
    SelectionModifiedObserver->SetClientData(this);

    PartSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionModifiedObserver
    );
    FieldSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionModifiedObserver
    );
    lagrangianFieldSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionModifiedObserver
    );
    LagrangianFieldSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionModifiedObserver
    );

    // Setup the animation-property-modified callback to reset the play-mode
    // back to to SNAP_TO_TIMESTEPS whenever ParaView mysteriously changes it
    // back to SEQUENCE
    vtkSMProxy* animationSceneProxy = getActiveAnimationSceneProxy();
    if (!animationSceneProxy) return;

    AnimationPropertyModifiedObserver = vtkCallbackCommand::New();
    AnimationPropertyModifiedObserver->SetCallback
    (
        &vtkPVFoamReader::AnimationPropertyModifiedCallback
    );
    AnimationPropertyModifiedObserver->SetClientData(this);

    animationSceneProxy->AddObserver
    (
        vtkCommand::PropertyModifiedEvent,
        this->AnimationPropertyModifiedObserver
    );

    // Explicitly call the animation-property-modified callback so that the
    // play-mode is correctly set at the start
    AnimationPropertyModifiedCallback(animationSceneProxy, 0, this, nullptr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

vtkPVFoamReader::~vtkPVFoamReader()
{
    vtkDebugMacro(<<"Deconstructor");

    if (foamData_)
    {
        // Remove patch names
        updatePatchNamesView(false);
        delete foamData_;
    }

    if (FileName)
    {
        delete [] FileName;
    }

    PartSelection->RemoveObserver(SelectionModifiedObserver);
    FieldSelection->RemoveObserver(SelectionModifiedObserver);
    lagrangianFieldSelection->RemoveObserver(SelectionModifiedObserver);
    LagrangianFieldSelection->RemoveObserver(SelectionModifiedObserver);

    SelectionModifiedObserver->Delete();

    PartSelection->Delete();
    FieldSelection->Delete();
    lagrangianFieldSelection->Delete();
    LagrangianFieldSelection->Delete();

    vtkSMProxy* animationSceneProxy = getActiveAnimationSceneProxy();
    if (!animationSceneProxy) return;

    animationSceneProxy->RemoveObserver(AnimationPropertyModifiedObserver);

    AnimationPropertyModifiedObserver->Delete();
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

    if (!FileName)
    {
        vtkErrorMacro("FileName has to be specified!");
        return 0;
    }

    if (Foam::vtkPVFoam::debug)
    {
        int nInfo = outputVector->GetNumberOfInformationObjects();

        cout<< "\nRequestInformation with " << nInfo << " item(s)\n";
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

    if (!requestTimeSteps(outputVector))
    {
        vtkErrorMacro("could not find valid OpenFOAM mesh");

        // delete foamData and flag it as fatal error
        delete foamData_;
        foamData_ = nullptr;
        return 0;
    }
    else
    {
        return 1;
    }
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

    if (First)
    {
        First = false;
        requestTimeSteps(outputVector);
    }

    int nInfo = outputVector->GetNumberOfInformationObjects();

    if (Foam::vtkPVFoam::debug)
    {
        cout<< "\nRequestData with " << nInfo << " item(s)\n";
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
        cout<< "Update output with "
            << output->GetNumberOfBlocks() << " blocks\n";
    }

    foamData_->Update(output, output, output);

    updatePatchNamesView(ShowPatchNames);

    // Do any cleanup on the OpenFOAM side
    foamData_->CleanUp();

    return 1;
}


void vtkPVFoamReader::SetRefresh()
{
    Modified();

    pqApplicationCore* appCore = pqApplicationCore::instance();
    if (!appCore) return;

    appCore->render();
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


void vtkPVFoamReader::PrintSelf(ostream& os, vtkIndent indent)
{
    vtkDebugMacro(<<"PrintSelf");

    this->Superclass::PrintSelf(os,indent);
    os  << indent << "File name: "
        << (this->FileName ? this->FileName : "(none)") << "\n";

    foamData_->PrintSelf(os, indent);

    os  << indent << "Time step: " << this->GetTimeStep() << endl;
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
// Field selection list control

vtkDataArraySelection* vtkPVFoamReader::GetFieldSelection()
{
    vtkDebugMacro(<<"GetFieldSelection");
    return FieldSelection;
}


int vtkPVFoamReader::GetNumberOfFieldArrays()
{
    vtkDebugMacro(<<"GetNumberOfFieldArrays");
    return FieldSelection->GetNumberOfArrays();
}


const char* vtkPVFoamReader::GetFieldArrayName(int index)
{
    vtkDebugMacro(<<"GetFieldArrayName");
    return FieldSelection->GetArrayName(index);
}


int vtkPVFoamReader::GetFieldArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetFieldArrayStatus");
    return FieldSelection->ArrayIsEnabled(name);
}


void vtkPVFoamReader::SetFieldArrayStatus(const char* name, int status)
{
    vtkDebugMacro(<<"SetFieldArrayStatus");
    if (status)
    {
        FieldSelection->EnableArray(name);
    }
    else
    {
        FieldSelection->DisableArray(name);
    }
}


// ----------------------------------------------------------------------
// lagrangianField selection list control

vtkDataArraySelection* vtkPVFoamReader::GetlagrangianFieldSelection()
{
    vtkDebugMacro(<<"GetlagrangianFieldSelection");
    return lagrangianFieldSelection;
}


int vtkPVFoamReader::GetNumberOflagrangianFieldArrays()
{
    vtkDebugMacro(<<"GetNumberOflagrangianFieldArrays");
    return lagrangianFieldSelection->GetNumberOfArrays();
}


const char* vtkPVFoamReader::GetlagrangianFieldArrayName(int index)
{
    vtkDebugMacro(<<"GetlagrangianFieldArrayName");
    return lagrangianFieldSelection->GetArrayName(index);
}


int vtkPVFoamReader::GetlagrangianFieldArrayStatus(const char* name)
{
    vtkDebugMacro(<<"GetlagrangianFieldArrayStatus");
    return lagrangianFieldSelection->ArrayIsEnabled(name);
}


void vtkPVFoamReader::SetlagrangianFieldArrayStatus
(
    const char* name,
    int status
)
{
    vtkDebugMacro(<<"SetlagrangianFieldArrayStatus");
    if (status)
    {
        lagrangianFieldSelection->EnableArray(name);
    }
    else
    {
        lagrangianFieldSelection->DisableArray(name);
    }
}


// ----------------------------------------------------------------------
// LagrangianField selection list control

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
    static_cast<vtkPVFoamReader*>(clientdata)->Modified();
}


void vtkPVFoamReader::AnimationPropertyModifiedCallback
(
    vtkObject* object,
    unsigned long,
    void* clientdata,
    void*
)
{
    vtkSMPropertyHelper
    (
        static_cast<vtkSMAnimationSceneProxy*>(object),
        "PlayMode"
    ).Set(vtkCompositeAnimationPlayer::SNAP_TO_TIMESTEPS);
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
