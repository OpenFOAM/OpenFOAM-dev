/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "blockMeshCartesianConfiguration.H"
#include "dictionary.H"
#include "polyPatch.H"
#include "wallPolyPatch.H"
#include "blockMeshFunctions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::List<Foam::word> Foam::blockMeshCartesianConfiguration::patches =
    {"xMin", "xMax", "yMin", "yMax", "zMin", "zMax"};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::blockMeshCartesianConfiguration::calcBlockMeshDict
(
    const bool& boundsOpt
)
{
    Info<< "Surface bounding box is " << bb_ << endl;

    // Round the bounding box if it is not specified with '-bounds' option
    const scalar roundFactor = roundingScale(bb_.minDim());
    if (!boundsOpt)
    {
        roundBoundingBox(bb_, roundFactor);
    }

    // Set nCells with the lowest number of cells within the range 10-20
    if (nCells_ == Vector<label>::zero)
    {
        nCells_ = Vector<label>(bb_.span()/roundFactor);

        if (!isEven(nCells_) && cmptMin(nCells_) > 20 && !boundsOpt)
        {
            roundBoundingBox(bb_, 2*roundFactor);
            nCells_ = Vector<label>(bb_.span()/roundFactor);
        }

        if (cmptMin(nCells_) < 10)
        {
            nCells_ *= 2;
        }
        else if (cmptMin(nCells_) > 20)
        {
            nCells_ /= 2;
        }
    }

    Info<< "Bounding box is now " << bb_ << endl;

    // Scale nCells_ by refine factor
    nCells_ *= refineFactor_;

    Info<< "Using background mesh nCells " << nCells_ << endl;
}


void Foam::blockMeshCartesianConfiguration::writeBackgroundMesh()
{
    dictionary dict("backgroundMesh");

    dict.add("xMin", bb_.min().x(), true);
    dict.add("xMax", bb_.max().x(), true);
    dict.add("yMin", bb_.min().y(), true);
    dict.add("yMax", bb_.max().y(), true);
    dict.add("zMin", bb_.min().z(), true);
    dict.add("zMax", bb_.max().z(), true);
    dict.add("xCells", nCells_.x(), true);
    dict.add("yCells", nCells_.y(), true);
    dict.add("zCells", nCells_.z(), true);

    os_ << dict.name().c_str()
        << dict << nl
        << "convertToMeters 1;" << nl
        << endl;
}


void Foam::blockMeshCartesianConfiguration::writeDefaultPatch()
{
    Pair<word> defaultPatch;

    word opt = "defaultPatch";
    if (patchOpts_.found(opt))
    {
        defaultPatch = readPatchOption(opt);
    }
    else
    {
        defaultPatch = {"background", "internal"};
    }

    beginDict(os_, "defaultPatch");

    os_ << indent << "name " << defaultPatch.first() << ";" << nl
        << indent << "type " << defaultPatch.second() << ";" << endl;

    endDict(os_);

    Info<< "\nAdding defaultPatch '" << defaultPatch.first()
        << "' of type '" << defaultPatch.second() << "'" << endl;
}


void Foam::blockMeshCartesianConfiguration::writePatch
(
    const word& name,
    const word& type,
    const string& face
)
{
    os_ << indent << name
        << " { type " << type
        << "; faces ( " << face.c_str()
        << " ); }" << endl;

    Info<< "Adding patch '" << name
        << "' of type '" << type << "'" << endl;
}


void Foam::blockMeshCartesianConfiguration::writeBoundary()
{
    // Enable boundary section if a patch option or clearBoundary is selected
    bool enableBoundary = clearBoundary_;
    forAll(patches, i)
    {
        if (enableBoundary)
        {
            break;
        }

        enableBoundary = patchOpts_.found(patches[i] + "Patch");
    }

    if (!enableBoundary)
    {
        os_ << "// delete \"-disabled\" to enable boundary settings" << endl;

        Info<< "\nNote: The boundary list in blockMeshDict is disabled" << nl
           << "To enable, open the file and edit line number "
           << os_.lineNumber() << nl << endl;
    }

    beginList
    (
        os_,
        enableBoundary ? "boundary" : "boundary-disabled"
    );

    const List<word> faces
    {
        "(0 3 7 4)",
        "(1 5 6 2)",
        "(0 4 5 1)",
        "(3 2 6 7)",
        "(0 1 2 3)",
        "(4 7 6 5)"
    };

    forAll(patches, i)
    {
        const bool optFound(patchOpts_.found(patches[i] + "Patch"));

        // Do not write patch entry if clearBoundary option is selected
        // and the respective patch option is not selected
        if (clearBoundary_ && !optFound)
        {
            continue;
        }

        Pair<word> patch(patches[i], "patch");

        if (optFound)
        {
            patch = readPatchOption(patch.first() + "Patch");
        }

        writePatch(patch.first(), patch.second(), faces[i]);
    }

    endList(os_);
}


void Foam::blockMeshCartesianConfiguration::writeVertices()
{
    beginList(os_, "vertices");

    writeVertex("xMin", "yMin", "zMin");
    writeVertex("xMax", "yMin", "zMin");
    writeVertex("xMax", "yMax", "zMin");
    writeVertex("xMin", "yMax", "zMin");
    writeVertex("xMin", "yMin", "zMax");
    writeVertex("xMax", "yMin", "zMax");
    writeVertex("xMax", "yMax", "zMax");
    writeVertex("xMin", "yMax", "zMax");

    endList(os_);
}


void Foam::blockMeshCartesianConfiguration::writeBlocks()
{
    beginList(os_, "blocks");

    os_ << indent << "hex (0 1 2 3 4 5 6 7)" << nl
        << indent << "(" << incrIndent << nl
        << indent << "$!backgroundMesh/xCells" << nl
        << indent << "$!backgroundMesh/yCells" << nl
        << indent << "$!backgroundMesh/zCells" << decrIndent << nl
        << indent << ")" << nl
        << indent << "simpleGrading (1 1 1)" << endl;

    endList(os_);
}


void Foam::blockMeshCartesianConfiguration::writeEdges()
{
    beginList(os_, "edges");
    endList(os_);
}


void Foam::blockMeshCartesianConfiguration::writeMergePatchPairs()
{
    beginList(os_, "mergePatchPairs");
    endList(os_, false);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockMeshCartesianConfiguration::blockMeshCartesianConfiguration
(
    const fileName& name,
    const fileName& dir,
    const Time& time,
    const meshingSurfaceList& surfaces,
    const bool& boundsOpt,
    const Vector<label>& nCells,
    const label refineFactor,
    const HashTable<Pair<word>>& patchOpts,
    const bool clearBoundary
)
:
    blockMeshConfigurationBase(name, dir, time, surfaces, patchOpts),
    nCells_(nCells),
    refineFactor_(refineFactor),
    clearBoundary_(clearBoundary)
{
    calcBlockMeshDict(boundsOpt);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blockMeshCartesianConfiguration::~blockMeshCartesianConfiguration()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::blockMeshCartesianConfiguration::write()
{
    dict_.writeHeader(os_, word("dictionary"));

    writeBackgroundMesh();
    writeDefaultPatch();
    writeBoundary();
    writeVertices();
    writeBlocks();
    writeEdges();
    writeMergePatchPairs();

    dict_.writeEndDivider(os_);
}


// ************************************************************************* //
