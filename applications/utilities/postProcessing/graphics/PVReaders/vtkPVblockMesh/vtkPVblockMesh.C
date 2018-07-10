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

#include "vtkPVblockMesh.H"
#include "vtkPVblockMeshReader.h"

// OpenFOAM includes
#include "blockMesh.H"
#include "Time.H"
#include "patchZones.H"
#include "OStringStream.H"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vtkPVblockMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtkPVblockMesh::resetCounters()
{
    // Reset mesh part ids and sizes
    arrayRangeBlocks_.reset();
    arrayRangeEdges_.reset();
    arrayRangeCorners_.reset();
}


void Foam::vtkPVblockMesh::updateInfoBlocks
(
    vtkDataArraySelection* arraySelection
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPVblockMesh::updateInfoBlocks"
            << " [meshPtr=" << (meshPtr_ ? "set" : "nullptr") << "]" << endl;
    }

    arrayRangeBlocks_.reset( arraySelection->GetNumberOfArrays() );

    const blockMesh& blkMesh = *meshPtr_;

    const int nBlocks = blkMesh.size();
    for (int blockI = 0; blockI < nBlocks; ++blockI)
    {
        const blockDescriptor& blockDef = blkMesh[blockI];

        // Display either blockI as a number or with its name
        // (looked up from blockMeshDict)
        OStringStream os;
        blockDescriptor::write(os, blockI, blkMesh.meshDict());
        word partName(os.str());

        // append the (optional) zone name
        if (!blockDef.zoneName().empty())
        {
            partName += " - " + blockDef.zoneName();
        }

        // Add blockId and zoneName to GUI list
        arraySelection->AddArray(partName.c_str());
    }

    arrayRangeBlocks_ += nBlocks;

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(arraySelection);

        Info<< "<end> Foam::vtkPVblockMesh::updateInfoBlocks" << endl;
    }
}


void Foam::vtkPVblockMesh::updateInfoEdges
(
    vtkDataArraySelection* arraySelection
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPVblockMesh::updateInfoEdges"
            << " [meshPtr=" << (meshPtr_ ? "set" : "nullptr") << "]" << endl;
    }

    arrayRangeEdges_.reset( arraySelection->GetNumberOfArrays() );

    const blockMesh& blkMesh = *meshPtr_;
    const blockEdgeList& edges = blkMesh.edges();

    const int nEdges = edges.size();
    forAll(edges, edgeI)
    {
        OStringStream ostr;
        blockVertex::write(ostr, edges[edgeI].start(), blkMesh.meshDict());
        ostr<< ":";
        blockVertex::write(ostr, edges[edgeI].end(), blkMesh.meshDict());
        ostr << " - " << edges[edgeI].type();

        // Add "beg:end - type" to GUI list
        arraySelection->AddArray(ostr.str().c_str());
    }

    arrayRangeEdges_ += nEdges;

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(arraySelection);

        Info<< "<end> Foam::vtkPVblockMesh::updateInfoEdges" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkPVblockMesh::vtkPVblockMesh
(
    const char* const FileName,
    vtkPVblockMeshReader* reader
)
:
    reader_(reader),
    dbPtr_(nullptr),
    meshPtr_(nullptr),
    meshRegion_(polyMesh::defaultRegion),
    meshDir_(polyMesh::meshSubDir),
    arrayRangeBlocks_("block"),
    arrayRangeEdges_("edges"),
    arrayRangeCorners_("corners")
{
    if (debug)
    {
        Info<< "Foam::vtkPVblockMesh::vtkPVblockMesh - "
            << FileName << endl;
    }

    // avoid argList and get rootPath/caseName directly from the file
    fileName fullCasePath(fileName(FileName).path());

    if (!isDir(fullCasePath))
    {
        return;
    }
    if (fullCasePath == ".")
    {
        fullCasePath = cwd();
    }

    // Set the case as an environment variable - some BCs might use this
    if (fullCasePath.name().find("processor", 0) == 0)
    {
        const fileName globalCase = fullCasePath.path();

        setEnv("FOAM_CASE", globalCase, true);
        setEnv("FOAM_CASENAME", globalCase.name(), true);
    }
    else
    {
        setEnv("FOAM_CASE", fullCasePath, true);
        setEnv("FOAM_CASENAME", fullCasePath.name(), true);
    }

    // look for 'case{region}.OpenFOAM'
    // could be stringent and insist the prefix match the directory name...
    // Note: cannot use fileName::name() due to the embedded '{}'
    string caseName(fileName(FileName).lessExt());
    string::size_type beg = caseName.find_last_of("/{");
    string::size_type end = caseName.find('}', beg);

    if
    (
        beg != string::npos && caseName[beg] == '{'
     && end != string::npos && end == caseName.size()-1
    )
    {
        meshRegion_ = caseName.substr(beg+1, end-beg-1);

        // some safety
        if (meshRegion_.empty())
        {
            meshRegion_ = polyMesh::defaultRegion;
        }

        if (meshRegion_ != polyMesh::defaultRegion)
        {
            meshDir_ = meshRegion_/polyMesh::meshSubDir;
        }
    }

    if (debug)
    {
        Info<< "fullCasePath=" << fullCasePath << nl
            << "FOAM_CASE=" << getEnv("FOAM_CASE") << nl
            << "FOAM_CASENAME=" << getEnv("FOAM_CASENAME") << endl;
    }

    // Create time object
    dbPtr_.reset
    (
        new Time
        (
            Time::controlDictName,
            fileName(fullCasePath.path()),
            fileName(fullCasePath.name())
        )
    );

    dbPtr_().functionObjects().off();

    updateInfo();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkPVblockMesh::~vtkPVblockMesh()
{
    if (debug)
    {
        Info<< "<end> Foam::vtkPVblockMesh::~vtkPVblockMesh" << endl;
    }

    // Hmm. pointNumberTextActors are not getting removed
    //
    forAll(pointNumberTextActorsPtrs_, pointi)
    {
        pointNumberTextActorsPtrs_[pointi]->Delete();
    }
    pointNumberTextActorsPtrs_.clear();

    delete meshPtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPVblockMesh::updateInfo()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPVblockMesh::updateInfo"
            << " [meshPtr=" << (meshPtr_ ? "set" : "nullptr") << "] " << endl;
    }

    resetCounters();

    vtkDataArraySelection* blockSelection = reader_->GetBlockSelection();
    vtkDataArraySelection* edgeSelection = reader_->GetCurvedEdgesSelection();

    // enable 'internalMesh' on the first call
    // or preserve the enabled selections
    stringList enabledParts;
    stringList enabledEdges;
    bool firstTime = false;
    if (!blockSelection->GetNumberOfArrays() && !meshPtr_)
    {
        firstTime = true;
    }
    else
    {
        enabledParts = getSelectedArrayEntries(blockSelection);
        enabledEdges = getSelectedArrayEntries(edgeSelection);
    }

    // Clear current mesh parts list
    blockSelection->RemoveAllArrays();
    edgeSelection->RemoveAllArrays();

    // need a blockMesh
    updateFoamMesh();

    // Update mesh parts list
    updateInfoBlocks( blockSelection );

    // Update curved edges list
    updateInfoEdges( edgeSelection );

    // restore the enabled selections
    if (!firstTime)
    {
        setSelectedArrayEntries(blockSelection, enabledParts);
        setSelectedArrayEntries(edgeSelection, enabledEdges);
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPVblockMesh::updateInfo" << endl;
    }
}


void Foam::vtkPVblockMesh::updateFoamMesh()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPVblockMesh::updateFoamMesh" << endl;
    }

    // Check to see if the OpenFOAM mesh has been created
    if (!meshPtr_)
    {
        if (debug)
        {
            Info<< "Creating blockMesh at time=" << dbPtr_().timeName()
                << endl;
        }

        // Set path for the blockMeshDict
        const word dictName("blockMeshDict");
        fileName dictPath(dbPtr_().system()/dictName);

        // Check if dictionary is present in the constant directory
        if
        (
            exists
            (
                dbPtr_().path()/dbPtr_().constant()
               /polyMesh::meshSubDir/dictName
            )
        )
        {
            dictPath = dbPtr_().constant()/polyMesh::meshSubDir/dictName;
        }

        // Store dictionary since is used as database inside blockMesh class
        // for names of vertices and blocks
        IOdictionary* meshDictPtr = new IOdictionary
        (
            IOobject
            (
                dictPath,
                dbPtr_(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                true
            )
        );
        meshDictPtr->store();

        meshPtr_ = new blockMesh(*meshDictPtr, meshRegion_);
    }


    if (debug)
    {
        Info<< "<end> Foam::vtkPVblockMesh::updateFoamMesh" << endl;
    }
}


void Foam::vtkPVblockMesh::Update
(
    vtkMultiBlockDataSet* output
)
{
    reader_->UpdateProgress(0.1);

    // Set up mesh parts selection(s)
    updateBoolListStatus(blockStatus_, reader_->GetBlockSelection());

    // Set up curved edges selection(s)
    updateBoolListStatus(edgeStatus_, reader_->GetCurvedEdgesSelection());

    reader_->UpdateProgress(0.2);

    // Update the OpenFOAM mesh
    updateFoamMesh();
    reader_->UpdateProgress(0.5);

    // Convert mesh element
    int blockNo = 0;

    convertMeshCorners(output, blockNo);
    convertMeshBlocks(output, blockNo);
    convertMeshEdges(output, blockNo);

    reader_->UpdateProgress(0.8);

}


void Foam::vtkPVblockMesh::CleanUp()
{
    reader_->UpdateProgress(1.0);
}


void Foam::vtkPVblockMesh::renderPointNumbers
(
    vtkRenderer* renderer,
    const bool show
)
{
    // always remove old actors first

    forAll(pointNumberTextActorsPtrs_, pointi)
    {
        renderer->RemoveViewProp(pointNumberTextActorsPtrs_[pointi]);
        pointNumberTextActorsPtrs_[pointi]->Delete();
    }
    pointNumberTextActorsPtrs_.clear();

    if (show && meshPtr_)
    {
        const blockMesh& blkMesh = *meshPtr_;
        const pointField& cornerPts = blkMesh.vertices();
        const scalar scaleFactor = blkMesh.scaleFactor();

        pointNumberTextActorsPtrs_.setSize(cornerPts.size());
        forAll(cornerPts, pointi)
        {
            vtkTextActor* txt = vtkTextActor::New();

            // Display either pointi as a number or with its name
            // (looked up from blockMeshDict)
            {
                OStringStream os;
                blockVertex::write(os, pointi, blkMesh.meshDict());
                txt->SetInput(os.str().c_str());
            }

            // Set text properties
            vtkTextProperty* tprop = txt->GetTextProperty();
            tprop->SetFontFamilyToArial();
            tprop->BoldOn();
            tprop->ShadowOff();
            tprop->SetLineSpacing(1.0);
            tprop->SetFontSize(14);
            tprop->SetColor(1.0, 0.0, 1.0);
            tprop->SetJustificationToCentered();

            // Set text to use 3-D world co-ordinates
            txt->GetPositionCoordinate()->SetCoordinateSystemToWorld();

            txt->GetPositionCoordinate()->SetValue
            (
                cornerPts[pointi].x()*scaleFactor,
                cornerPts[pointi].y()*scaleFactor,
                cornerPts[pointi].z()*scaleFactor
            );

            // Add text to each renderer
            renderer->AddViewProp(txt);

            // Maintain a list of text labels added so that they can be
            // removed later
            pointNumberTextActorsPtrs_[pointi] = txt;
        }
    }
}



void Foam::vtkPVblockMesh::PrintSelf(ostream& os, vtkIndent indent) const
{
#if 0
    os  << indent << "Number of nodes: "
        << (meshPtr_ ? meshPtr_->nPoints() : 0) << "\n";

    os  << indent << "Number of cells: "
        << (meshPtr_ ? meshPtr_->nCells() : 0) << "\n";

    os  << indent << "Number of available time steps: "
        << (dbPtr_.valid() ? dbPtr_().times().size() : 0) << endl;
#endif
}

// ************************************************************************* //
