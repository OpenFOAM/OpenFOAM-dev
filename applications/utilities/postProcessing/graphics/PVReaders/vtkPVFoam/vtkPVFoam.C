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

#include "vtkPVFoam.H"
#include "vtkPVFoamReader.h"

// OpenFOAM includes
#include "domainDecomposition.H"
#include "fvFieldReconstructor.H"
#include "pointFieldReconstructor.H"
#include "lagrangianFieldReconstructor.H"
#include "LagrangianFieldReconstructor.H"
#include "fvMeshStitcher.H"
#include "patchZones.H"
#include "fileOperation.H"
#include "IFstream.H"
#include "OSspecific.H"
#include "etcFiles.H"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vtkPVFoam, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtkPVFoam::clearReconstructors()
{
    fvReconstructorPtr_.clear();
    pointReconstructorPtr_.clear();
    lagrangianReconstructors_.clear();
    LagrangianMeshes_.clear();
    LagrangianReconstructors_.clear();
}


void Foam::vtkPVFoam::clearFoamMesh()
{
    if
    (
        !reader_->GetCacheMesh()
     || (
            procMeshesPtr_.valid()
         && reader_->GetDecomposedCase() != procMeshesPtr_->haveProcs()
        )
    )
    {
        clearReconstructors();

        procMeshesPtr_.clear();
    }
}


void Foam::vtkPVFoam::resetCounters()
{
    arrayRangeVolume_.reset();
    arrayRangePatches_.reset();
    arrayRangelagrangian_.reset();
    arrayRangeLagrangian_.reset();
    arrayRangeCellZones_.reset();
    arrayRangeFaceZones_.reset();
    arrayRangePointZones_.reset();
    arrayRangeCellSets_.reset();
    arrayRangeFaceSets_.reset();
    arrayRangePointSets_.reset();
}


void Foam::vtkPVFoam::reduceMemory()
{
    forAll(regionPolyDecomp_, i)
    {
        regionPolyDecomp_[i].clear();
    }

    forAll(zonePolyDecomp_, i)
    {
        zonePolyDecomp_[i].clear();
    }

    forAll(setPolyDecomp_, i)
    {
        setPolyDecomp_[i].clear();
    }

    clearFoamMesh();
}


int Foam::vtkPVFoam::setTime(int nRequest, const double requestTimes[])
{
    DebugInFunction;

    const Time& runTime =
        reader_->GetDecomposedCase()
      ? procDbsPtr_->proc0Time()
      : procDbsPtr_->completeTime();

    // Get times list. Flush first to force refresh.
    fileHandler().flush();
    const instantList times = runTime.times();

    // Find the nearest index to the selected time
    int nearestIndex = timeIndex_;
    for (int requestI = 0; requestI < nRequest; ++requestI)
    {
        int index = Time::findClosestTimeIndex(times, requestTimes[requestI]);
        if (index >= 0 && index != timeIndex_)
        {
            nearestIndex = index;
            break;
        }
    }

    // Clip the index to zero
    if (nearestIndex < 0) nearestIndex = 0;

    // If the time has changed...
    if (timeIndex_ != nearestIndex)
    {
        // Set the time
        timeIndex_ = nearestIndex;
        procDbsPtr_->setTime(times[nearestIndex], nearestIndex);

        // Clear the mesh if necessary
        clearFoamMesh();

        // Update the mesh
        fvMesh::readUpdateState stat = fvMesh::TOPO_PATCH_CHANGE;
        if (procMeshesPtr_.valid())
        {
            if (reader_->GetDecomposedCase())
            {
                stat = procMeshesPtr_->readUpdateReconstruct(true);
            }
            else
            {
                stat = procMeshesPtr_->readUpdateComplete();
            }
        }

        if (stat > fvMesh::POINTS_MOVED)
        {
            clearReconstructors();
        }
    }

    // Re-stitch if necessary
    if (procMeshesPtr_.valid())
    {
        procMeshesPtr_->completeMesh().stitcher().reconnect
        (
            reader_->GetInterpolateVolFields()
        );
    }

    return nearestIndex;
}


void Foam::vtkPVFoam::topoChangePartsStatus()
{
    vtkDataArraySelection* selection = reader_->GetPartSelection();

    const label nElem = selection->GetNumberOfArrays();

    // Clear the part statuses
    if (partStatus_.size() != nElem)
    {
        partStatus_.setSize(nElem);
        partStatus_ = false;
    }

    // Clear the part datasets. Note that this is not optimal as it means we
    // are not re-using existing data sets.
    partDataset_.setSize(nElem);
    partDataset_ = -1;

    // Read the selected mesh parts (zones, patches ...) and add to list
    forAll(partStatus_, partId)
    {
        const int setting = selection->GetArraySetting(partId);

        if (partStatus_[partId] != setting)
        {
            partStatus_[partId] = setting;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkPVFoam::vtkPVFoam
(
    const char* const FileNameCStr,
    vtkPVFoamReader* reader
)
:
    reader_(reader),
    procDbsPtr_(nullptr),
    procMeshesPtr_(nullptr),
    meshRegion_(polyMesh::defaultRegion),
    meshDir_(polyMesh::meshSubDir),
    timeIndex_(-1),
    arrayRangeVolume_("unzoned"),
    arrayRangePatches_("patches"),
    arrayRangelagrangian_("lagrangian"),
    arrayRangeLagrangian_("Lagrangian"),
    arrayRangeCellZones_("cellZone"),
    arrayRangeFaceZones_("faceZone"),
    arrayRangePointZones_("pointZone"),
    arrayRangeCellSets_("cellSet"),
    arrayRangeFaceSets_("faceSet"),
    arrayRangePointSets_("pointSet")
{
    DebugInFunction
        << "fileName=" << FileNameCStr << endl;

    fileName FileName(FileNameCStr);

    // Avoid argList and get rootPath/caseName directly from the file
    fileName fullCasePath(FileName.path());

    if (!isDir(fullCasePath)) return;

    if (fullCasePath == ".")
    {
        fullCasePath = cwd();
    }

    if (fullCasePath.name().find("processors", 0) == 0)
    {
        // FileName e.g. "cavity/processors256/processor1.OpenFOAM
        // Remove the processors section so it goes into processorDDD
        // checking below.
        fullCasePath = fullCasePath.path()/fileName(FileName.name()).lessExt();
    }

    if (fullCasePath.name().find("processor", 0) == 0)
    {
        // Give filehandler opportunity to analyse number of processors
        (void)fileHandler().filePath(fullCasePath);

        const fileName globalCase = fullCasePath.path();

        setEnv("FOAM_CASE", globalCase, true);
        setEnv("FOAM_CASENAME", globalCase.name(), true);
    }
    else
    {
        setEnv("FOAM_CASE", fullCasePath, true);
        setEnv("FOAM_CASENAME", fullCasePath.name(), true);
    }

    // Parse the mesh and region names from 'case(mesh){region}' in FileName
    string caseName(FileName.name(true));

    string::size_type beg = caseName.find_last_of('(');
    string::size_type end = caseName.find(')', beg);

    if (beg != string::npos && end != string::npos)
    {
        meshMesh_ = caseName.substr(beg+1, end-beg-1);
        meshPath_ = "meshes"/meshMesh_;
        meshDir_ = meshPath_/polyMesh::meshSubDir;
    }

    beg = caseName.find_last_of('{');
    end = caseName.find('}', beg);

    if (beg != string::npos && end != string::npos)
    {
        meshRegion_ = caseName.substr(beg+1, end-beg-1);
        meshDir_ = meshPath_/meshRegion_/polyMesh::meshSubDir;
    }

    DebugInfo
        << "    fullCasePath=" << fullCasePath << nl
        << "    FOAM_CASE=" << getEnv("FOAM_CASE") << nl
        << "    FOAM_CASENAME=" << getEnv("FOAM_CASENAME") << nl
        << "    mesh=" << meshMesh_ << nl
        << "    region=" << meshRegion_ << endl;

    // Pre-load any libraries
    dlLibraryTable dlTable;
    string libsString(getEnv("FOAM_LIBS"));
    if (!libsString.empty())
    {
        IStringStream is(libsString);
        fileNameList libNames(is);
        forAll(libNames, i)
        {
            dlTable.open(libNames[i]);
        }
    }

    // Create time object
    procDbsPtr_.reset
    (
        new processorRunTimes
        (
            Time::controlDictName,
            fileName(fullCasePath.path()),
            fileName(fullCasePath.name()),
            false,
            processorRunTimes::nProcsFrom::fileHandler
        )
    );

    // Read the configuration
    fileNameList configDictFiles = findEtcFiles("paraFoam", false);
    forAllReverse(configDictFiles, cdfi)
    {
        configDict_.merge(dictionary(IFstream(configDictFiles[cdfi])()));
    }

    updateInfo();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkPVFoam::~vtkPVFoam()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPVFoam::updateInfo()
{
    DebugInFunction;

    resetCounters();

    vtkDataArraySelection* partSelection = reader_->GetPartSelection();

    // Determine whether or not this is the first update
    const bool first =
        !partSelection->GetNumberOfArrays() && !procMeshesPtr_.valid();

    // Enable 'internalMesh' on the first call, otherwise or preserve the
    // previously enabled selections
    stringList enabledEntries;
    if (first)
    {
        enabledEntries.setSize(1);
        enabledEntries[0] = "internalMesh";
    }
    else
    {
        enabledEntries = getSelectedArrayEntries(partSelection, false);
    }

    // Clear current mesh parts list
    partSelection->RemoveAllArrays();

    // Update mesh parts list - add Lagrangian at the bottom
    updateInfoInternalMesh(partSelection);
    updateInfoPatches(partSelection, enabledEntries, first);
    updateInfoSets(partSelection);
    updateInfoZones(partSelection);
    updateInfolagrangian(partSelection);
    updateInfoLagrangian(partSelection);

    // Restore the enabled selections
    setSelectedArrayEntries(partSelection, enabledEntries);

    if (debug) getSelectedArrayEntries(partSelection);

    // Update fields
    updateInfoFields();
    updateInfolagrangianFields();
    updateInfoLagrangianFields();
}


void Foam::vtkPVFoam::updateFoamMesh()
{
    DebugInFunction;

    // Clear the mesh if necessary
    clearFoamMesh();

    // Create the OpenFOAM mesh if it does not yet exist
    if (!procMeshesPtr_.valid())
    {
        const bool haveMeshMesh = !meshMesh_.empty();
        const bool haveMeshRegion = meshRegion_ != polyMesh::defaultRegion;

        DebugInfo
            << "Creating OpenFOAM mesh"
            << (haveMeshMesh ? " for mesh " + meshMesh_ : "").c_str()
            << (haveMeshMesh && haveMeshRegion ? " and" : "")
            << (haveMeshRegion ? " for region " + meshRegion_ : "").c_str()
            << " at time=" << procDbsPtr_().completeTime().name() << endl;

        procMeshesPtr_.reset
        (
            new domainDecomposition
            (
                procDbsPtr_(),
                meshPath_,
                meshRegion_
            )
        );

        if (reader_->GetDecomposedCase())
        {
            procMeshesPtr_->readReconstruct(true);
        }
        else
        {
            procMeshesPtr_->readComplete();
        }
    }
    else
    {
        DebugInfo
            << "Using existing OpenFOAM mesh" << endl;
    }

    // Stitch if necessary
    if (procMeshesPtr_.valid())
    {
        procMeshesPtr_->completeMesh().stitcher().reconnect
        (
            reader_->GetInterpolateVolFields()
        );
    }
}


void Foam::vtkPVFoam::Update
(
    vtkMultiBlockDataSet* output,
    vtkMultiBlockDataSet* lagrangianOutput,
    vtkMultiBlockDataSet* LagrangianOutput
)
{
    DebugInFunction;

    reader_->UpdateProgress(0.1);

    // Set up mesh parts selection(s)
    topoChangePartsStatus();

    reader_->UpdateProgress(0.15);

    // Update the OpenFOAM mesh
    updateFoamMesh();
    reader_->UpdateProgress(0.4);

    // Convert meshes - start port0 at block=0
    int blockNo = 0;

    convertMeshVolume(output, blockNo);
    convertMeshPatches(output, blockNo);

    reader_->UpdateProgress(0.6);

    if (reader_->GetIncludeZones())
    {
        convertMeshCellZones(output, blockNo);
        convertMeshFaceZones(output, blockNo);
        convertMeshPointZones(output, blockNo);

        reader_->UpdateProgress(0.65);
    }

    if (reader_->GetIncludeSets())
    {
        convertMeshCellSets(output, blockNo);
        convertMeshFaceSets(output, blockNo);
        convertMeshPointSets(output, blockNo);

        reader_->UpdateProgress(0.7);
    }

    convertMeshlagrangian(lagrangianOutput, blockNo);
    convertMeshLagrangian(LagrangianOutput, blockNo);

    reader_->UpdateProgress(0.8);

    // Update fields
    convertFields(output);
    convertlagrangianFields(lagrangianOutput);
    convertLagrangianFields(LagrangianOutput);

    reader_->UpdateProgress(0.95);
}


void Foam::vtkPVFoam::CleanUp()
{
    // Reclaim some memory
    reduceMemory();
    reader_->UpdateProgress(1.0);
}


double* Foam::vtkPVFoam::findTimes(const bool first, int& nTimeSteps)
{
    // Read the available times in a given database
    auto findTimesForRunTime = [this](const Time& runTime)
    {
        // Get all the times for this database
        const instantList timeLst = runTime.times();

        // Find the first time for which this mesh appears to exist
        label timei = 0;
        for (; timei < timeLst.size(); ++timei)
        {
            if
            (
                typeIOobject<pointIOField>
                (
                    "points",
                    timeLst[timei].name(),
                    meshDir_,
                    runTime
                ).headerOk()
            )
            {
                break;
            }
        }

        label nTimes = timeLst.size() - timei;

        // Skip "constant" time whenever possible
        if (timei == 0 && nTimes > 1)
        {
            if (timeLst[timei].name() == Time::constant())
            {
                ++timei;
                --nTimes;
            }
        }

        // Skip "0/" time if requested and possible
        if (nTimes > 1 && reader_->GetSkipZeroTime())
        {
            if (mag(timeLst[timei].value()) < small)
            {
                ++timei;
                --nTimes;
            }
        }

        return instantList(SubList<instant>(timeLst, nTimes, timei));
    };

    // Get a list of available instants
    instantList times;
    if (procDbsPtr_.valid())
    {
        // Get times from complete and/or processor databases. Use both if this
        // is the first execution.
        const instantList completeTimes =
            first || !reader_->GetDecomposedCase()
          ? findTimesForRunTime(procDbsPtr_->completeTime())
          : instantList();

        const instantList procTimes =
            first || reader_->GetDecomposedCase()
          ? findTimesForRunTime(procDbsPtr_->proc0Time())
          : instantList();

        // Merge the lists of times
        times.resize(completeTimes.size() + procTimes.size());
        label completeTimei = 0, procTimei = 0, timei = 0;
        while
        (
            completeTimei < completeTimes.size()
         && procTimei < procTimes.size()
        )
        {
            const bool completeNext =
                completeTimes[completeTimei].value()
              < procTimes[procTimei].value();

            const bool procNext =
                completeTimes[completeTimei].value()
              > procTimes[procTimei].value();

            times[timei ++] =
                completeNext
              ? completeTimes[completeTimei]
              : procTimes[procTimei];

            if (!procNext) completeTimei ++;
            if (!completeNext) procTimei ++;
        }
        while (completeTimei < completeTimes.size())
        {
            times[timei ++] = completeTimes[completeTimei ++];
        }
        while (procTimei < procTimes.size())
        {
            times[timei ++] = procTimes[procTimei ++];
        }
        times.resize(timei);
    }

    // If we have some times, convert to a bare array for VTK and return
    if (times.size())
    {
        nTimeSteps = times.size();
        double* timeSteps = new double[times.size()];
        forAll(times, timei)
        {
            timeSteps[timei] = times[timei].value();
        }
        return timeSteps;
    }
    else
    {
        nTimeSteps = 0;
        return nullptr;
    }
}


void Foam::vtkPVFoam::renderPatchNames
(
    vtkRenderer* renderer,
    const bool show
)
{
    if (!procMeshesPtr_.valid()) return;

    // Always remove old actors first
    forAll(patchTextActorsPtrs_, patchi)
    {
        renderer->RemoveViewProp(patchTextActorsPtrs_[patchi]);
        patchTextActorsPtrs_[patchi]->Delete();
    }
    patchTextActorsPtrs_.clear();

    if (show)
    {
        // Get the display patches, strip off any suffix
        const wordHashSet selectedPatches = getSelected
        (
            reader_->GetPartSelection(),
            arrayRangePatches_
        );

        if (selectedPatches.empty()) return;

        const polyBoundaryMesh& pbMesh =
            procMeshesPtr_->completeMesh().boundaryMesh();

        // Find the total number of zones
        // Each zone will take the patch name
        // Number of zones per patch ... zero zones should be skipped
        labelList nZones(pbMesh.size(), 0);

        // Per global zone number the average face centre position
        List<DynamicList<point>> zoneCentre(pbMesh.size());

        // Loop through all patches to determine zones, and centre of each zone
        forAll(pbMesh, patchi)
        {
            const polyPatch& pp = pbMesh[patchi];

            // Only include the patch if it is selected
            if (!selectedPatches.found(pp.name())) continue;

            const labelListList& edgeFaces = pp.edgeFaces();
            const vectorField& n = pp.faceNormals();

            boolList featEdge(pp.nEdges(), false);

            forAll(edgeFaces, edgeI)
            {
                const labelList& eFaces = edgeFaces[edgeI];

                if (eFaces.size() == 1)
                {
                    // Note: could also do ones with > 2 faces but this gives
                    // too many zones for baffles
                    featEdge[edgeI] = true;
                }
                else if (mag(n[eFaces[0]] & n[eFaces[1]]) < 0.5)
                {
                    featEdge[edgeI] = true;
                }
            }

            // Do topological analysis of patch, find disconnected regions
            patchZones pZones(pp, featEdge);

            nZones[patchi] = pZones.nZones();

            labelList zoneNFaces(pZones.nZones(), 0);

            // Create storage for additional zone centres
            forAll(zoneNFaces, zoneI)
            {
                zoneCentre[patchi].append(Zero);
            }

            // Do averaging per individual zone
            forAll(pp, facei)
            {
                label zoneI = pZones[facei];
                zoneCentre[patchi][zoneI] += pp[facei].centre(pp.points());
                zoneNFaces[zoneI]++;
            }

            forAll(zoneCentre[patchi], zoneI)
            {
                zoneCentre[patchi][zoneI] /= zoneNFaces[zoneI];
            }
        }

        // Count number of zones we're actually going to display.
        // This is truncated to a max per patch

        const label MAXPATCHZONES = 20;

        label displayZoneI = 0;

        forAll(pbMesh, patchi)
        {
            displayZoneI += min(MAXPATCHZONES, nZones[patchi]);
        }

        // Set the size of the patch labels to max number of zones
        patchTextActorsPtrs_.setSize(displayZoneI);

        // Actor index
        displayZoneI = 0;

        forAll(pbMesh, patchi)
        {
            const polyPatch& pp = pbMesh[patchi];

            label globalZoneI = 0;

            // Only selected patches will have a non-zero number of zones
            label nDisplayZones = min(MAXPATCHZONES, nZones[patchi]);
            label increment = 1;
            if (nZones[patchi] >= MAXPATCHZONES)
            {
                increment = nZones[patchi]/MAXPATCHZONES;
            }

            for (label i = 0; i < nDisplayZones; i++)
            {
                vtkTextActor* txt = vtkTextActor::New();

                txt->SetInput(pp.name().c_str());

                // Set text properties
                vtkTextProperty* tprop = txt->GetTextProperty();
                tprop->SetFontFamilyToArial();
                tprop->BoldOff();
                tprop->ShadowOff();
                tprop->SetLineSpacing(1.0);
                tprop->SetFontSize(12);
                tprop->SetColor(1.0, 0.0, 0.0);
                tprop->SetJustificationToCentered();

                // Set text to use 3-D world co-ordinates
                txt->GetPositionCoordinate()->SetCoordinateSystemToWorld();

                txt->GetPositionCoordinate()->SetValue
                (
                    zoneCentre[patchi][globalZoneI].x(),
                    zoneCentre[patchi][globalZoneI].y(),
                    zoneCentre[patchi][globalZoneI].z()
                );

                // Add text to each renderer
                renderer->AddViewProp(txt);

                // Maintain a list of text labels added so that they can be
                // removed later
                patchTextActorsPtrs_[displayZoneI] = txt;

                globalZoneI += increment;
                displayZoneI++;
            }
        }

        // Resize the patch names list to the actual number of patch names added
        patchTextActorsPtrs_.setSize(displayZoneI);
    }
}


void Foam::vtkPVFoam::PrintSelf(ostream& os, vtkIndent indent) const
{
    os  << indent << "Number of nodes: " << (procMeshesPtr_.valid()
         ? procMeshesPtr_->completeMesh().nPoints() : 0) << "\n";

    os  << indent << "Number of cells: " << (procMeshesPtr_.valid()
         ? procMeshesPtr_->completeMesh().nCells() : 0) << "\n";

    const label nTimeSteps =
        procDbsPtr_.empty()
      ? 0
      : reader_->GetDecomposedCase()
      ? procDbsPtr_->proc0Time().times().size()
      : procDbsPtr_->completeTime().times().size();

    os  << indent << "Number of available time steps: " << nTimeSteps << "\n";

    os  << indent << "mesh: " << meshMesh_ << "\n";

    os  << indent << "mesh region: " << meshRegion_ << "\n";
}


// ************************************************************************* //
