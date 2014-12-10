/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "directions.H"
#include "polyMesh.H"
#include "twoDPointCorrector.H"
#include "directionInfo.H"
#include "MeshWave.H"
#include "OFstream.H"
#include "meshTools.H"
#include "hexMatcher.H"
#include "Switch.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::directions::directionType,
        3
    >::names[] =
    {
        "tan1",
        "tan2",
        "normal"
    };
}

const Foam::NamedEnum<Foam::directions::directionType, 3>
    Foam::directions::directionTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// For debugging
void Foam::directions::writeOBJ(Ostream& os, const point& pt)
{
    os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
}

// For debugging
void Foam::directions::writeOBJ
(
    Ostream& os,
    const point& pt0,
    const point& pt1,
    label& vertI
)
{
    writeOBJ(os, pt0);
    writeOBJ(os, pt1);

    os << "l " << vertI + 1 << ' ' << vertI + 2 << endl;

    vertI += 2;
}


// Dump to file.
void Foam::directions::writeOBJ
(
    const fileName& fName,
    const primitiveMesh& mesh,
    const vectorField& dirs
)
{
    Pout<< "Writing cell info to " << fName << " as vectors at the cellCentres"
        << endl << endl;

    OFstream xDirStream(fName);

    label vertI = 0;

    forAll(dirs, cellI)
    {
        const point& ctr = mesh.cellCentres()[cellI];

        // Calculate local length scale
        scalar minDist = GREAT;

        const labelList& nbrs = mesh.cellCells()[cellI];

        forAll(nbrs, nbrI)
        {
            minDist = min(minDist, mag(mesh.cellCentres()[nbrs[nbrI]] - ctr));
        }

        scalar scale = 0.5*minDist;

        writeOBJ(xDirStream, ctr, ctr + scale*dirs[cellI], vertI);
    }
}


void Foam::directions::check2D
(
    const twoDPointCorrector* correct2DPtr,
    const vector& vec
)
{
    if (correct2DPtr)
    {
        if (mag(correct2DPtr->planeNormal() & vec) > 1e-6)
        {
            FatalErrorIn("check2D") << "Specified vector " << vec
                << "is not normal to plane defined in dynamicMeshDict."
                << endl
                << "Either make case 3D or adjust vector."
                << exit(FatalError);
        }
    }
}


// Get direction on all cells
Foam::vectorField Foam::directions::propagateDirection
(
    const polyMesh& mesh,
    const bool useTopo,
    const polyPatch& pp,
    const vectorField& ppField,
    const vector& defaultDir
)
{
    // Seed all faces on patch
    labelList changedFaces(pp.size());
    List<directionInfo> changedFacesInfo(pp.size());

    if (useTopo)
    {
        forAll(pp, patchFaceI)
        {
            label meshFaceI = pp.start() + patchFaceI;

            label cellI = mesh.faceOwner()[meshFaceI];

            if (!hexMatcher().isA(mesh, cellI))
            {
                FatalErrorIn("propagateDirection")
                    << "useHexTopology specified but cell " << cellI
                    << " on face " << patchFaceI << " of patch " << pp.name()
                    << " is not a hex" << exit(FatalError);
            }

            const vector& cutDir = ppField[patchFaceI];

            // Get edge(bundle) on cell most in direction of cutdir
            label edgeI = meshTools::cutDirToEdge(mesh, cellI, cutDir);

            // Convert edge into index on face
            label faceIndex =
                directionInfo::edgeToFaceIndex
                (
                    mesh,
                    cellI,
                    meshFaceI,
                    edgeI
                );

            // Set initial face and direction
            changedFaces[patchFaceI] = meshFaceI;
            changedFacesInfo[patchFaceI] =
                directionInfo
                (
                    faceIndex,
                    cutDir
                );
        }
    }
    else
    {
        forAll(pp, patchFaceI)
        {
            changedFaces[patchFaceI] = pp.start() + patchFaceI;
            changedFacesInfo[patchFaceI] =
                directionInfo
                (
                    -2,         // Geometric information only
                    ppField[patchFaceI]
                );
        }
    }

    MeshWave<directionInfo> directionCalc
    (
        mesh,
        changedFaces,
        changedFacesInfo,
        mesh.globalData().nTotalCells()+1
    );

    const List<directionInfo>& cellInfo = directionCalc.allCellInfo();

    vectorField dirField(cellInfo.size());

    label nUnset = 0;
    label nGeom = 0;
    label nTopo = 0;

    forAll(cellInfo, cellI)
    {
        label index = cellInfo[cellI].index();

        if (index == -3)
        {
            // Never visited
            WarningIn("propagateDirection")
                << "Cell " << cellI << " never visited to determine "
                << "local coordinate system" << endl
                << "Using direction " << defaultDir << " instead" << endl;

            dirField[cellI] = defaultDir;

            nUnset++;
        }
        else if (index == -2)
        {
            // Geometric direction
            dirField[cellI] = cellInfo[cellI].n();

            nGeom++;
        }
        else if (index == -1)
        {
            FatalErrorIn("propagateDirection")
                << "Illegal index " << index << endl
                << "Value is only allowed on faces" << abort(FatalError);
        }
        else
        {
            // Topological edge cut. Convert into average cut direction.
            dirField[cellI] = meshTools::edgeToCutDir(mesh, cellI, index);

            nTopo++;
        }
    }

    Pout<< "Calculated local coords for " << defaultDir
        << endl
        << "    Geometric cut cells   : " << nGeom << endl
        << "    Topological cut cells : " << nTopo << endl
        << "    Unset cells           : " << nUnset << endl
        << endl;

    return dirField;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::directions::directions
(
    const polyMesh& mesh,
    const dictionary& dict,
    const twoDPointCorrector* correct2DPtr
)
:
    List<vectorField>(wordList(dict.lookup("directions")).size())
{
    const wordList wantedDirs(dict.lookup("directions"));

    bool wantNormal = false;
    bool wantTan1 = false;
    bool wantTan2 = false;

    forAll(wantedDirs, i)
    {
        directionType wantedDir = directionTypeNames_[wantedDirs[i]];

        if (wantedDir == NORMAL)
        {
            wantNormal = true;
        }
        else if (wantedDir == TAN1)
        {
            wantTan1 = true;
        }
        else if (wantedDir == TAN2)
        {
            wantTan2 = true;
        }
    }


    label nDirs = 0;

    const word coordSystem(dict.lookup("coordinateSystem"));

    if (coordSystem == "global")
    {
        const dictionary& globalDict = dict.subDict("globalCoeffs");

        vector tan1(globalDict.lookup("tan1"));
        check2D(correct2DPtr, tan1);

        vector tan2(globalDict.lookup("tan2"));
        check2D(correct2DPtr, tan2);

        vector normal = tan1 ^ tan2;
        normal /= mag(normal);

        Pout<< "Global Coordinate system:" << endl
            << "     normal : " << normal << endl
            << "     tan1   : " << tan1 << endl
            << "     tan2   : " << tan2
            << endl << endl;

        if (wantNormal)
        {
            operator[](nDirs++) = vectorField(1, normal);
        }
        if (wantTan1)
        {
            operator[](nDirs++) = vectorField(1, tan1);
        }
        if (wantTan2)
        {
            operator[](nDirs++) = vectorField(1, tan2);
        }
    }
    else if (coordSystem == "patchLocal")
    {
        const dictionary& patchDict = dict.subDict("patchLocalCoeffs");

        const word patchName(patchDict.lookup("patch"));

        const label patchI = mesh.boundaryMesh().findPatchID(patchName);

        if (patchI == -1)
        {
            FatalErrorIn
            (
                "directions::directions(const polyMesh&, const dictionary&,"
                "const twoDPointCorrector*)"
            )   << "Cannot find patch "
                << patchName
                << exit(FatalError);
        }

        // Take zeroth face on patch
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        vector tan1(patchDict.lookup("tan1"));

        const vector& n0 = pp.faceNormals()[0];

        if (correct2DPtr)
        {
            tan1 = correct2DPtr->planeNormal() ^ n0;

            WarningIn
            (
                "directions::directions(const polyMesh&, const dictionary&,"
                "const twoDPointCorrector*)"
            )   << "Discarding user specified tan1 since 2D case." << endl
                << "Recalculated tan1 from face normal and planeNormal as "
                << tan1 << endl << endl;
        }

        Switch useTopo(dict.lookup("useHexTopology"));

        vectorField normalDirs;
        vectorField tan1Dirs;

        if (wantNormal || wantTan2)
        {
            normalDirs =
                propagateDirection
                (
                    mesh,
                    useTopo,
                    pp,
                    pp.faceNormals(),
                    n0
                );

            if (wantNormal)
            {
                this->operator[](nDirs++) = normalDirs;
            }
        }

        if (wantTan1 || wantTan2)
        {
            tan1Dirs =
                propagateDirection
                (
                    mesh,
                    useTopo,
                    pp,
                    vectorField(pp.size(), tan1),
                    tan1
                );


            if (wantTan1)
            {
                this->operator[](nDirs++) = tan1Dirs;
            }
        }
        if (wantTan2)
        {
            tmp<vectorField> tan2Dirs = normalDirs ^ tan1Dirs;

            this->operator[](nDirs++) = tan2Dirs;
        }
    }
    else
    {
        FatalErrorIn
        (
            "directions::directions(const polyMesh&, const dictionary&,"
            "const twoDPointCorrector*)"
        )   << "Unknown coordinate system "
            << coordSystem << endl
            << "Known types are global and patchLocal"
            << exit(FatalError);
    }
}


// ************************************************************************* //
