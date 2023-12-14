/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "faceCoupleInfo.H"
#include "polyMesh.H"
#include "matchPoints.H"
#include "indirectPrimitivePatch.H"
#include "meshTools.H"
#include "treeDataFace.H"
#include "indexedOctree.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceCoupleInfo, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceCoupleInfo::writeOBJ
(
    const fileName& fName,
    const edgeList& edges,
    const pointField& points
)
{
    OFstream str(fName);

    labelList pointMap(points.size(), -1);

    label newPointi = 0;

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        forAll(e, eI)
        {
            label pointi = e[eI];

            if (pointMap[pointi] == -1)
            {
                pointMap[pointi] = newPointi++;

                meshTools::writeOBJ(str, points[pointi]);
            }
        }
    }

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        str<< "l " << pointMap[e[0]]+1 << ' ' << pointMap[e[1]]+1 << nl;
    }
}


void Foam::faceCoupleInfo::writeOBJ
(
    const fileName& fName,
    const pointField& points0,
    const pointField& points1
)
{
    Pout<< "Writing connections as edges to " << fName << endl;

    OFstream str(fName);

    label vertI = 0;

    forAll(points0, i)
    {
        meshTools::writeOBJ(str, points0[i]);
        vertI++;
        meshTools::writeOBJ(str, points1[i]);
        vertI++;
        str << "l " << vertI-1 << ' ' << vertI << nl;
    }
}


void Foam::faceCoupleInfo::writePointsFaces() const
{
    const indirectPrimitivePatch& m = masterPatch();
    const indirectPrimitivePatch& s = slavePatch();

    // Patches
    {
        OFstream str("masterPatch.obj");
        Pout<< "Writing masterPatch to " << str.name() << endl;
        meshTools::writeOBJ(str, m.localFaces(), m.localPoints());
    }
    {
        OFstream str("slavePatch.obj");
        Pout<< "Writing slavePatch to " << str.name() << endl;
        meshTools::writeOBJ(str, s.localFaces(), s.localPoints());
    }

    // Point connectivity
    {
        Pout<< "Writing masterToSlavePoints to masterToSlavePoints.obj" << endl;

        const labelListList coupleToMasterPoints(this->coupleToMasterPoints());
        const labelListList coupleToSlavePoints(this->coupleToSlavePoints());
        pointField coupleMasterPoints(coupleToMasterPoints.size());
        pointField coupleSlavePoints(coupleToSlavePoints.size());
        forAll(coupleToMasterPoints, couplePointi)
        {
            const label masterPointi = coupleToMasterPoints[couplePointi][0];
            const label slavePointi = coupleToSlavePoints[couplePointi][0];
            coupleMasterPoints[couplePointi] = m.localPoints()[masterPointi];
            coupleSlavePoints[couplePointi] = m.localPoints()[slavePointi];
        }

        writeOBJ
        (
            "masterToSlavePoints.obj",
            coupleMasterPoints,
            coupleSlavePoints
        );
    }

    // Face connectivity
    {
        Pout<< "Writing masterToSlaveFaces to masterToSlaveFaces.obj" << endl;

        writeOBJ
        (
            "masterToSlaveFaces.obj",
            calcFaceCentres<IndirectList>(m, m.points(), 0, m.size()),
            calcFaceCentres<IndirectList>(s, s.points(), 0, s.size())
        );
    }

    Pout<< endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceCoupleInfo::faceCoupleInfo
(
    const polyMesh& masterMesh,
    const labelList& masterAddressing,
    const polyMesh& slaveMesh,
    const labelList& slaveAddressing
)
:
    masterPatchPtr_
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>(masterMesh.faces(), masterAddressing),
            masterMesh.points()
        )
    ),
    slavePatchPtr_
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>(slaveMesh.faces(), slaveAddressing),
            slaveMesh.points()
        )
    ),
    nCouplePoints_(0),
    masterToCouplePoints_(),
    slaveToCouplePoints_()
{
    if (masterAddressing.size() != slaveAddressing.size())
    {
        FatalErrorInFunction
            << "Number of master and slave faces differ." << endl
            << "master:" << masterAddressing.size()
            << " slave:" << slaveAddressing.size()
            << abort(FatalError);
    }

    if
    (
        masterAddressing.size()
     && min(masterAddressing) < masterMesh.nInternalFaces()
    )
    {
        FatalErrorInFunction
            << "Supplied internal face on master mesh to couple." << nl
            << "Faces to be coupled have to be boundary faces."
            << abort(FatalError);
    }

    if
    (
        slaveAddressing.size()
     && min(slaveAddressing) < slaveMesh.nInternalFaces()
    )
    {
        FatalErrorInFunction
            << "Supplied internal face on slave mesh to couple." << nl
            << "Faces to be coupled have to be boundary faces."
            << abort(FatalError);
    }

    // Initialise point addressing
    nCouplePoints_ = 0;
    masterToCouplePoints_ = labelList(masterPatch().nPoints(), -1);
    slaveToCouplePoints_ = labelList(slavePatch().nPoints(), -1);

    // Slave points number around the master faces in the reverse direction and
    // have the same first point
    forAll(masterPatch(), coupleFacei)
    {
        const face& masterF = masterPatch().localFaces()[coupleFacei];
        const face& slaveF = slavePatch().localFaces()[coupleFacei];

        label slaveFp = 0;

        forAll(masterF, masterFp)
        {
            const label masterPointi = masterF[masterFp];
            const label slavePointi = slaveF[slaveFp];

            // If this point already has a coupled point index then use it,
            // else create a new one
            label couplePointi = -1;
            if (masterToCouplePoints_[masterPointi] != -1)
            {
                couplePointi = masterToCouplePoints_[masterPointi];
            }
            else if (slaveToCouplePoints_[slavePointi] != -1)
            {
                couplePointi = slaveToCouplePoints_[slavePointi];
            }
            else
            {
                couplePointi = nCouplePoints_ ++;
            }

            masterToCouplePoints_[masterPointi] = couplePointi;
            slaveToCouplePoints_[slavePointi] = couplePointi;

            slaveFp = slaveF.rcIndex(slaveFp);
        }
    }

    if (debug)
    {
        writePointsFaces();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceCoupleInfo::~faceCoupleInfo()
{}


// ************************************************************************* //
