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

#include "blockMeshCylindricalConfiguration.H"
#include "dictionary.H"
#include "polyPatch.H"
#include "wallPolyPatch.H"
#include "unitConversion.H"
#include "blockMeshFunctions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::List<Foam::word> Foam::blockMeshCylindricalConfiguration::patches =
    {"zMin", "zMax"};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::blockMeshCylindricalConfiguration::bbInflate
(
    boundBox& bb,
    const vector& s
)
{
    bb.min() = cmptMultiply(s, bb.min());
    bb.max() = cmptMultiply(s, bb.max());
}

bool Foam::blockMeshCylindricalConfiguration::isBoundBoxOnZaxis()
{
    const vector bbMidNorm = bb_.midpoint()/bb_.mag();
    return mag(bbMidNorm.x()) < rootSmall && mag(bbMidNorm.y()) < rootSmall;
}

void Foam::blockMeshCylindricalConfiguration::calcBlockMeshDict()
{
    // Set nCells as a vector: (boxCells radialCells zCells)
    label boxCells = nCells_.x();
    label radialCells = nCells_.y();

    if (nCells_ == Vector<label>::zero)
    {
        boxCells = 5;
        radialCells = 15;

        nCells_ = Vector<label>
        (
            boxCells,
            radialCells,
            2.0*radialCells*bb_.span().z()/bb_.span().x()
        );
    }

    // Size the bounding box
    const scalar roundFactor = roundingScale(0.01*bb_.minDim());
    const scalar expansion = 1.0/(Foam::cos(degToRad(45.0/nCells_.x())));

    const vector scaling(expansion, expansion, 1);

    // Inflate the bounding box in radial direction
    bbInflate(bb_, scaling);
    bbInflate(rzbb_, scaling);

    // Round the bounding box
    roundBoundingBox(bb_, roundFactor);
    roundBoundingBox(rzbb_, roundFactor);

    radBox_ = ceil(0.3*rzbb_.max().x()/roundFactor)*roundFactor;

    nCells_ *= refineFactor_;
}


void Foam::blockMeshCylindricalConfiguration::writeBackgroundMesh()
{
    const scalar radOut = bb_.max().x();
    const scalar radIn = rzbb_.max().x();
    const scalar boxToRadOut = radOut - radBox_;
    const scalar boxToRadIn = radIn - radBox_;

    const label inCells = ceil(boxToRadIn*nCells_.y()/boxToRadOut);
    const label outCells = nCells_.y() - inCells;

    beginDict(os_, "backgroundMesh");

    os_ << indent << "radOut     " << radOut << ";" << endl;
    os_ << indent << "radIn      " << radIn << ";" << endl;
    os_ << indent << "radBox     " << radBox_ << ";" << nl << endl;
    os_ << indent << "zMin       " << bb_.min().z() << ";" << endl;
    os_ << indent << "zMax       " << bb_.max().z() << ";" << nl << endl;
    os_ << indent << "boxCells   " << nCells_.x() << ";" << endl;
    os_ << indent << "inCells    " << inCells << ";" << endl;
    os_ << indent << "outCells   " << outCells << ";" << endl;
    os_ << indent << "zCells     " << nCells_.z() << ";" << nl << endl;
    os_ << indent << "outGrading 2.0;" << nl << endl;
    os_ << indent << "radOutN    #neg $radOut;" << endl;
    os_ << indent << "radInN     #neg $radIn;" << endl;
    os_ << indent << "radBoxN    #neg $radBox;" << endl;

    endDict(os_);

    os_ << "convertToMeters 1;" << nl << endl;
}


void Foam::blockMeshCylindricalConfiguration::writeDefaultPatch()
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


void Foam::blockMeshCylindricalConfiguration::writePatch
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


void Foam::blockMeshCylindricalConfiguration::writeBoundary()
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

        beginDict(os_, patch.first());
        os_ << indent << "type " << patch.second() << ";" << endl;
        beginList(os_, "faces");

        switch (i)
        {
            case 0:
            {
                os_ << indent << "(0 1 2 3)" << nl
                    << indent << "(0 4 5 1)" << nl
                    << indent << "(1 5 6 2)" << nl
                    << indent << "(2 6 7 3)" << nl
                    << indent << "(3 7 4 0)" << nl
                    << indent << "(4 8 9 5)" << nl
                    << indent << "(5 9 10 6)" << nl
                    << indent << "(6 10 11 7)" << nl
                    << indent << "(7 11 8 4)" << endl;
                break;
            }

            case 1:
            {
                os_ << indent << "(12 13 14 15)" << nl
                    << indent << "(12 16 17 13)" << nl
                    << indent << "(13 17 18 14)" << nl
                    << indent << "(14 18 19 15)" << nl
                    << indent << "(15 19 16 12)" << nl
                    << indent << "(16 20 21 17)" << nl
                    << indent << "(17 21 22 18)" << nl
                    << indent << "(18 22 23 19)" << nl
                    << indent << "(19 23 20 16)" << endl;
                break;
            }
        }

        endList(os_, false);
        endDict(os_, i != 1);
    }

    endList(os_);
}


void Foam::blockMeshCylindricalConfiguration::writeGeometry()
{
    beginDict(os_, "geometry");

    List<word> geometries {"rotatingZone", "outer"};
    List<word> dims {"radIn", "radOut"};

    const scalar zMin = roundDown(bb_.min().z(), 10);
    const scalar zMax = roundUp(bb_.max().z(), 10);

    forAll(geometries, i)
    {
        beginDict(os_, geometries[i]);

        os_ << indent << "type searchableCylinder;" << nl
            << indent << "point1 (0 0 " << zMin << ");" << nl
            << indent << "point2 (0 0 " << zMax << ");" << nl
            << indent << "radius $!backgroundMesh/" << dims[i] << ";" << endl;

        endDict(os_);
    }

    endDict(os_);
}


void Foam::blockMeshCylindricalConfiguration::writeProjectedVertex
(
    const word& x,
    const word& y,
    const word& z,
    const word& surface
)
{
    os_ << indent << "project" << endl;
    writeVertex(x, y, z);
    os_ << indent << "(" << surface << ")" << nl << endl;
}


void Foam::blockMeshCylindricalConfiguration::writeVertices()
{
    beginList(os_, "vertices");

    forAll(patches, i)
    {
        const word& dir = patches[i];

        writeVertex("radBoxN", "radBoxN", dir);
        writeVertex("radBox", "radBoxN", dir);
        writeVertex("radBox", "radBox", dir);
        writeVertex("radBoxN", "radBox", dir);

        os_ << endl;

        writeProjectedVertex("radInN", "radInN", dir, "rotatingZone");
        writeProjectedVertex("radIn", "radInN", dir, "rotatingZone");
        writeProjectedVertex("radIn", "radIn", dir, "rotatingZone");
        writeProjectedVertex("radInN", "radIn", dir, "rotatingZone");

        writeProjectedVertex("radOutN", "radOutN", dir, "outer");
        writeProjectedVertex("radOut", "radOutN", dir, "outer");
        writeProjectedVertex("radOut", "radOut", dir, "outer");
        writeProjectedVertex("radOutN", "radOut", dir, "outer");
    }

    endList(os_);
}


void Foam::blockMeshCylindricalConfiguration::writeBlocks()
{
    os_ << "boxMesh" << endl;
    writeVertex("boxCells", "boxCells", "zCells");
    os_ << "simpleGrading (1 1 1);" << nl << endl;

    os_ << "inMesh" << endl;
    writeVertex("boxCells", "inCells", "zCells");
    os_ << "simpleGrading (1 1 1);" << nl << endl;

    os_ << "outMesh" << endl;
    writeVertex("boxCells", "outCells", "zCells");
    os_ << "simpleGrading (1 $!backgroundMesh/outGrading 1);" << nl << endl;

    beginList(os_, "blocks");

    os_ << indent << "hex (0 1 2 3 12 13 14 15) $boxMesh" << nl << nl
        << indent << "hex (1 0 4 5 13 12 16 17) $inMesh" << nl
        << indent << "hex (0 3 7 4 12 15 19 16) $inMesh" << nl
        << indent << "hex (2 1 5 6 14 13 17 18) $inMesh" << nl
        << indent << "hex (3 2 6 7 15 14 18 19) $inMesh" << nl << nl
        << indent << "hex (5 4  8  9 17 16 20 21) $outMesh" << nl
        << indent << "hex (4 7 11  8 16 19 23 20) $outMesh" << nl
        << indent << "hex (6 5  9 10 18 17 21 22) $outMesh" << nl
        << indent << "hex (7 6 10 11 19 18 22 23) $outMesh" << endl;

    endList(os_);
}


void Foam::blockMeshCylindricalConfiguration::writeEdges()
{
    beginList(os_, "edges");

    os_ << indent << "project  4  5  (rotatingZone)" << nl
        << indent << "project  5  6  (rotatingZone)" << nl
        << indent << "project  6  7  (rotatingZone)" << nl
        << indent << "project  7  4  (rotatingZone)" << nl
        << indent << "project 16 17  (rotatingZone)" << nl
        << indent << "project 17 18  (rotatingZone)" << nl
        << indent << "project 18 19  (rotatingZone)" << nl
        << indent << "project 19 16  (rotatingZone)" << nl << nl
        << indent << "project  8  9  (outer)" << nl
        << indent << "project  9 10  (outer)" << nl
        << indent << "project 10 11  (outer)" << nl
        << indent << "project 11  8  (outer)" << nl
        << indent << "project 20 21  (outer)" << nl
        << indent << "project 21 22  (outer)" << nl
        << indent << "project 22 23  (outer)" << nl
        << indent << "project 23 20  (outer)" << endl;

    endList(os_);
}


void Foam::blockMeshCylindricalConfiguration::writeMergePatchPairs()
{
    os_ << "mergePatchPairs" << nl
        << "(" << nl
        << ");" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockMeshCylindricalConfiguration::blockMeshCylindricalConfiguration
(
    const fileName& name,
    const fileName& dir,
    const Time& time,
    const meshingSurfaceList& surfaces,
    const Vector<label>& nCells,
    const label refineFactor,
    const HashTable<Pair<word>>& patchOpts,
    const bool clearBoundary
)
:
    blockMeshConfigurationBase(name, dir, time, surfaces, patchOpts),
    rzbb_(surfaces.rzbb()),
    nCells_(nCells),
    refineFactor_(refineFactor),
    clearBoundary_(clearBoundary)
{
    if (!isBoundBoxOnZaxis())
    {
        FatalErrorInFunction
            << "Attempting to create a cylindrical background mesh"
            << nl << "but the geometry bounds are not aligned with the z-axis."
            << exit(FatalError);
    }

    if (rzbb_.volume() == 0)
    {
        WarningInFunction
            << "Creating a cylindrical background mesh without a "
            << "rotatingZone specified by the '-rotatingZones' option."
            << nl <<endl;

        // Create the intermediate interface at 40% of domain size if no
        // rotating zone is specified
        scalar factor = 0.4;

        rzbb_.min() = factor*bb_.min();
        rzbb_.max() = factor*bb_.max();
    }

    calcBlockMeshDict();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blockMeshCylindricalConfiguration::~blockMeshCylindricalConfiguration()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::blockMeshCylindricalConfiguration::write()
{
    dict_.writeHeader(os_, word("dictionary"));

    writeBackgroundMesh();
    writeDefaultPatch();
    writeBoundary();
    writeGeometry();
    writeVertices();
    writeBlocks();
    writeEdges();
    writeMergePatchPairs();

    dict_.writeEndDivider(os_);
}


// ************************************************************************* //
