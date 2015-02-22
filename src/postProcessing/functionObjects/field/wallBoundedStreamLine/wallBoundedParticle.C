/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "wallBoundedParticle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const std::size_t Foam::wallBoundedParticle::sizeofFields_
(
    sizeof(wallBoundedParticle) - sizeof(particle)
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::edge Foam::wallBoundedParticle::currentEdge() const
{
    if ((meshEdgeStart_ != -1) == (diagEdge_ != -1))
    {
        FatalErrorIn("wallBoundedParticle::currentEdge() const")
            << "Particle:"
            << info()
            << "cannot both be on a mesh edge and a face-diagonal edge."
            << " meshEdgeStart_:" << meshEdgeStart_
            << " diagEdge_:" << diagEdge_
            << abort(FatalError);
    }

    const Foam::face& f = mesh_.faces()[tetFace()];

    if (meshEdgeStart_ != -1)
    {
        return edge(f[meshEdgeStart_], f.nextLabel(meshEdgeStart_));
    }
    else
    {
        label faceBasePtI = mesh_.tetBasePtIs()[tetFace()];
        label diagPtI = (faceBasePtI+diagEdge_)%f.size();

        return edge(f[faceBasePtI], f[diagPtI]);
    }
}


void Foam::wallBoundedParticle::crossEdgeConnectedFace
(
    const edge& meshEdge
)
{
    // Update tetFace, tetPt
    particle::crossEdgeConnectedFace(cell(), tetFace(), tetPt(), meshEdge);

    // Update face to be same as tracking one
    face() = tetFace();

    // And adapt meshEdgeStart_.
    const Foam::face& f = mesh_.faces()[tetFace()];
    label fp = findIndex(f, meshEdge[0]);

    if (f.nextLabel(fp) == meshEdge[1])
    {
        meshEdgeStart_ = fp;
    }
    else
    {
        label fpMin1 = f.rcIndex(fp);

        if (f[fpMin1] == meshEdge[1])
        {
            meshEdgeStart_ = fpMin1;
        }
        else
        {
            FatalErrorIn
            (
                "wallBoundedParticle::crossEdgeConnectedFace"
                "(const edge&)"
            )   << "Problem :"
                << " particle:"
                << info()
                << "face:" << tetFace()
                << " verts:" << f
                << " meshEdge:" << meshEdge
                << abort(FatalError);
        }
    }

    diagEdge_ = -1;

    // Check that still on same mesh edge
    const edge eNew(f[meshEdgeStart_], f.nextLabel(meshEdgeStart_));
    if (eNew != meshEdge)
    {
        FatalErrorIn
        (
            "wallBoundedParticle::crossEdgeConnectedFace"
            "(const edge&)"
        )   << "Problem" << abort(FatalError);
    }
}


void Foam::wallBoundedParticle::crossDiagonalEdge()
{
    if (diagEdge_ == -1)
    {
        FatalErrorIn("wallBoundedParticle::crossDiagonalEdge()")
            << "Particle:"
            << info()
            << "not on a diagonal edge" << abort(FatalError);
    }
    if (meshEdgeStart_ != -1)
    {
        FatalErrorIn("wallBoundedParticle::crossDiagonalEdge()")
            << "Particle:"
            << info()
            << "meshEdgeStart_:" << meshEdgeStart_ << abort(FatalError);
    }

    const Foam::face& f = mesh_.faces()[tetFace()];

    // tetPtI starts from 1, goes up to f.size()-2

    if (tetPt() == diagEdge_)
    {
        tetPt() = f.rcIndex(tetPt());
    }
    else
    {
        label nextTetPt = f.fcIndex(tetPt());
        if (diagEdge_ == nextTetPt)
        {
            tetPt() = nextTetPt;
        }
        else
        {
            FatalErrorIn("wallBoundedParticle::crossDiagonalEdge()")
                << "Particle:"
                << info()
                << "tetPt:" << tetPt()
                << " diagEdge:" << diagEdge_ << abort(FatalError);
        }
    }

    meshEdgeStart_ = -1;
}


Foam::scalar Foam::wallBoundedParticle::trackFaceTri
(
    const vector& endPosition,
    label& minEdgeI
)
{
    // Track p from position to endPosition
    const triFace tri(currentTetIndices().faceTriIs(mesh_));
    vector n = tri.normal(mesh_.points());
    n /= mag(n)+VSMALL;

    // Check which edge intersects the trajectory.
    // Project trajectory onto triangle
    minEdgeI = -1;
    scalar minS = 1;        // end position

    edge currentE(-1, -1);
    if (meshEdgeStart_ != -1 || diagEdge_ != -1)
    {
        currentE = currentEdge();
    }

    // Determine path along line position+s*d to see where intersections
    // are.

    forAll(tri, i)
    {
        label j = tri.fcIndex(i);

        const point& pt0 = mesh_.points()[tri[i]];
        const point& pt1 = mesh_.points()[tri[j]];

        if (edge(tri[i], tri[j]) == currentE)
        {
            // Do not check particle is on
            continue;
        }

        // Outwards pointing normal
        vector edgeNormal = (pt1-pt0)^n;

        edgeNormal /= mag(edgeNormal)+VSMALL;

        // Determine whether position and end point on either side of edge.
        scalar sEnd = (endPosition-pt0)&edgeNormal;
        if (sEnd >= 0)
        {
            // endPos is outside triangle. position() should always be
            // inside.
            scalar sStart = (position()-pt0)&edgeNormal;
            if (mag(sEnd-sStart) > VSMALL)
            {
                scalar s = sStart/(sStart-sEnd);

                if (s >= 0 && s < minS)
                {
                    minS = s;
                    minEdgeI = i;
                }
            }
        }
    }

    if (minEdgeI != -1)
    {
        position() += minS*(endPosition-position());
    }
    else
    {
        // Did not hit any edge so tracked to the end position. Set position
        // without any calculation to avoid truncation errors.
        position() = endPosition;
        minS = 1.0;
    }

    // Project position() back onto plane of triangle
    const point& triPt = mesh_.points()[tri[0]];
    position() -= ((position()-triPt)&n)*n;

    return minS;
}


bool Foam::wallBoundedParticle::isTriAlongTrack
(
    const point& endPosition
) const
{
    const triFace triVerts(currentTetIndices().faceTriIs(mesh_));
    const edge currentE = currentEdge();

    //if (debug)
    {
        if
        (
            currentE[0] == currentE[1]
         || findIndex(triVerts, currentE[0]) == -1
         || findIndex(triVerts, currentE[1]) == -1
        )
        {
            FatalErrorIn
            (
                "wallBoundedParticle::isTriAlongTrack"
                "(const point&)"
            )   << "Edge " << currentE << " not on triangle " << triVerts
                << info()
                << abort(FatalError);
        }
    }


    const vector dir = endPosition-position();

    // Get normal of currentE
    vector n = triVerts.normal(mesh_.points());
    n /= mag(n);

    forAll(triVerts, i)
    {
        label j = triVerts.fcIndex(i);
        const point& pt0 = mesh_.points()[triVerts[i]];
        const point& pt1 = mesh_.points()[triVerts[j]];

        if (edge(triVerts[i], triVerts[j]) == currentE)
        {
            vector edgeNormal = (pt1-pt0)^n;
            return (dir&edgeNormal) < 0;
        }
    }

    FatalErrorIn
    (
        "wallBoundedParticle::isTriAlongTrack"
        "(const point&)"
    )   << "Problem" << abort(FatalError);

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoundedParticle::wallBoundedParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label meshEdgeStart,
    const label diagEdge
)
:
    particle(mesh, position, cellI, tetFaceI, tetPtI),
    meshEdgeStart_(meshEdgeStart),
    diagEdge_(diagEdge)
{}


Foam::wallBoundedParticle::wallBoundedParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is  >> meshEdgeStart_ >> diagEdge_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&meshEdgeStart_),
                sizeof(meshEdgeStart_)
              + sizeof(diagEdge_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "wallBoundedParticle::wallBoundedParticle"
        "(const Cloud<wallBoundedParticle>&, Istream&, bool)"
    );
}


Foam::wallBoundedParticle::wallBoundedParticle
(
    const wallBoundedParticle& p
)
:
    particle(p),
    meshEdgeStart_(p.meshEdgeStart_),
    diagEdge_(p.diagEdge_)
{}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const wallBoundedParticle& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.meshEdgeStart_
            << token::SPACE << p.diagEdge_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.meshEdgeStart_),
            wallBoundedParticle::sizeofFields_
        );
    }

    return os;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<wallBoundedParticle>& ip
)
{
    const wallBoundedParticle& p = ip.t_;

    tetPointRef tpr(p.currentTet());

    os  << "    " << static_cast<const particle&>(p) << nl
        << "    tet:" << nl;
    os  << "    ";
    meshTools::writeOBJ(os, tpr.a());
    os  << "    ";
    meshTools::writeOBJ(os, tpr.b());
    os  << "    ";
    meshTools::writeOBJ(os, tpr.c());
    os  << "    ";
    meshTools::writeOBJ(os, tpr.d());
    os  << "    l 1 2" << nl
        << "    l 1 3" << nl
        << "    l 1 4" << nl
        << "    l 2 3" << nl
        << "    l 2 4" << nl
        << "    l 3 4" << nl;
    os  << "    ";
    meshTools::writeOBJ(os, p.position());

    return os;
}


// ************************************************************************* //
