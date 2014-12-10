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

#include "edgeStats.H"
#include "Time.H"
#include "polyMesh.H"
#include "Ostream.H"
#include "twoDPointCorrector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::edgeStats::edgeTol_ = 1e-3;



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::direction Foam::edgeStats::getNormalDir
(
    const twoDPointCorrector* correct2DPtr
) const
{
    direction dir = 3;

    if (correct2DPtr)
    {
        const vector& normal = correct2DPtr->planeNormal();

        if (mag(normal & vector(1, 0, 0)) > 1-edgeTol_)
        {
            dir = 0;
        }
        else if (mag(normal & vector(0, 1, 0)) > 1-edgeTol_)
        {
            dir = 1;
        }
        else if (mag(normal & vector(0, 0, 1)) > 1-edgeTol_)
        {
            dir = 2;
        }
    }
    return dir;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::edgeStats::edgeStats(const polyMesh& mesh)
:
    mesh_(mesh),
    normalDir_(3)
{
    IOobject motionObj
    (
        "motionProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    );

    if (motionObj.headerOk())
    {
        Info<< "Reading " << mesh.time().constant() / "motionProperties"
            << endl << endl;

        IOdictionary motionProperties(motionObj);

        Switch twoDMotion(motionProperties.lookup("twoDMotion"));

        if (twoDMotion)
        {
            Info<< "Correcting for 2D motion" << endl << endl;

            autoPtr<twoDPointCorrector> correct2DPtr
            (
                new twoDPointCorrector(mesh)
            );

            normalDir_ = getNormalDir(&correct2DPtr());
        }
    }
}


// Construct from components
Foam::edgeStats::edgeStats
(
    const polyMesh& mesh,
    const twoDPointCorrector* correct2DPtr
)
:
    mesh_(mesh),
    normalDir_(getNormalDir(correct2DPtr))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::edgeStats::minLen(Ostream& os) const
{
    label nX = 0;
    label nY = 0;
    label nZ = 0;

    scalar minX = GREAT;
    scalar maxX = -GREAT;
    vector x(1, 0, 0);

    scalar minY = GREAT;
    scalar maxY = -GREAT;
    vector y(0, 1, 0);

    scalar minZ = GREAT;
    scalar maxZ = -GREAT;
    vector z(0, 0, 1);

    scalar minOther = GREAT;
    scalar maxOther = -GREAT;

    const edgeList& edges = mesh_.edges();

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        vector eVec(e.vec(mesh_.points()));

        scalar eMag = mag(eVec);

        eVec /= eMag;

        if (mag(eVec & x) > 1-edgeTol_)
        {
            minX = min(minX, eMag);
            maxX = max(maxX, eMag);
            nX++;
        }
        else if (mag(eVec & y) > 1-edgeTol_)
        {
            minY = min(minY, eMag);
            maxY = max(maxY, eMag);
            nY++;
        }
        else if (mag(eVec & z) > 1-edgeTol_)
        {
            minZ = min(minZ, eMag);
            maxZ = max(maxZ, eMag);
            nZ++;
        }
        else
        {
            minOther = min(minOther, eMag);
            maxOther = max(maxOther, eMag);
        }
    }

    os  << "Mesh bounding box:" << boundBox(mesh_.points()) << nl << nl
        << "Mesh edge statistics:" << nl
        << "    x aligned :  number:" << nX << "\tminLen:" << minX
        << "\tmaxLen:" << maxX << nl
        << "    y aligned :  number:" << nY << "\tminLen:" << minY
        << "\tmaxLen:" << maxY << nl
        << "    z aligned :  number:" << nZ << "\tminLen:" << minZ
        << "\tmaxLen:" << maxZ << nl
        << "    other     :  number:" << mesh_.nEdges() - nX - nY - nZ
        << "\tminLen:" << minOther
        << "\tmaxLen:" << maxOther << nl << endl;

    if (normalDir_ == 0)
    {
        return min(minY, min(minZ, minOther));
    }
    else if (normalDir_ == 1)
    {
        return min(minX, min(minZ, minOther));
    }
    else if (normalDir_ == 2)
    {
        return min(minX, min(minY, minOther));
    }
    else
    {
        return min(minX, min(minY, min(minZ, minOther)));
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
