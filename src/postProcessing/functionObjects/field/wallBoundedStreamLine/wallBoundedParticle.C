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

#include "wallBoundedParticle.H"
#include "vectorFieldIOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
//    defineParticleTypeNameAndDebug(wallBoundedParticle, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//// Check position is inside tet
//void Foam::wallBoundedParticle::checkInside() const
//{
//    const tetIndices ti(currentTetIndices());
//    const tetPointRef tpr(ti.tet(mesh_));
//    if (!tpr.inside(position()))
//    {
//        FatalErrorIn("wallBoundedParticle::checkInside(..)")
//            << "Particle:" //<< static_cast<const particle&>(*this)
//            << info()
//            << "is not inside " << tpr
//            << abort(FatalError);
//    }
//}
//
//
//void Foam::wallBoundedParticle::checkOnEdge() const
//{
//    // Check that edge (as indicated by meshEdgeStart_, diagEdge_) is
//    // indeed one that contains the position.
//    const edge e = currentEdge();
//
//    linePointRef ln(e.line(mesh_.points()));
//
//    pointHit ph(ln.nearestDist(position()));
//
//    if (ph.distance() > 1e-6)
//    {
//        FatalErrorIn
//        (
//            "wallBoundedParticle::checkOnEdge()"
//        )   << "Problem :"
//            << " particle:" //<< static_cast<const particle&>(*this)
//            << info()
//            << "edge:" << e
//            << " at:" << ln
//            << " distance:" << ph.distance()
//            << abort(FatalError);
//    }
//}
//
//
//void Foam::wallBoundedParticle::checkOnTriangle(const point& p)
//const
//{
//    const triFace tri(currentTetIndices().faceTriIs(mesh_));
//    pointHit ph = tri.nearestPoint(p, mesh_.points());
//    if (ph.distance() > 1e-9)
//    {
//        FatalErrorIn
//        (
//            "wallBoundedParticle::checkOnTriangle(const point&)"
//        )   << "Problem :"
//            << " particle:" //<< static_cast<const particle&>(*this)
//            << info()
//            << "point:" << p
//            << " distance:" << ph.distance()
//            << abort(FatalError);
//    }
//}


// Construct the edge the particle is on (according to meshEdgeStart_,
// diagEdge_)
Foam::edge Foam::wallBoundedParticle::currentEdge() const
{
    if ((meshEdgeStart_ != -1) == (diagEdge_ != -1))
    {
        FatalErrorIn("wallBoundedParticle::currentEdge() const")
            << "Particle:" //<< static_cast<const particle&>(*this)
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
    //label oldFaceI = tetFace();

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
                << " particle:" //<< static_cast<const particle&>(*this)
                << info()
                << "face:" << tetFace()
                << " verts:" << f
                << " meshEdge:" << meshEdge
                << abort(FatalError);
        }
    }

    diagEdge_ = -1;

    //Pout<< "    crossed meshEdge "
    //    << meshEdge.line(mesh().points())
    //    << " from face:" << oldFaceI
    //    << " to face:" << tetFace() << endl;


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


    // Check that edge (as indicated by meshEdgeStart_) is indeed one that
    // contains the position.
    //checkOnEdge();
}


void Foam::wallBoundedParticle::crossDiagonalEdge()
{
    if (diagEdge_ == -1)
    {
        FatalErrorIn("wallBoundedParticle::crossDiagonalEdge()")
            << "Particle:" //<< static_cast<const particle&>(*this)
            << info()
            << "not on a diagonal edge" << abort(FatalError);
    }
    if (meshEdgeStart_ != -1)
    {
        FatalErrorIn("wallBoundedParticle::crossDiagonalEdge()")
            << "Particle:" //<< static_cast<const particle&>(*this)
            << info()
            << "meshEdgeStart_:" << meshEdgeStart_ << abort(FatalError);
    }

    //label oldTetPt = tetPt();

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
                << "Particle:" //<< static_cast<const particle&>(*this)
                << info()
                << "tetPt:" << tetPt()
                << " diagEdge:" << diagEdge_ << abort(FatalError);
        }
    }

    meshEdgeStart_ = -1;

    //Pout<< "    crossed diagonalEdge "
    //    << currentEdge().line(mesh().points())
    //    << " from tetPt:" << oldTetPt
    //    << " to tetPt:" << tetPt() << endl;
}


//- Track through a single triangle.
// Gets passed tet+triangle the particle is in. Updates position() but nothing
// else. Returns the triangle edge the particle is now on.
Foam::scalar Foam::wallBoundedParticle::trackFaceTri
(
    const vector& endPosition,
    label& minEdgeI
)
{
    // Track p from position to endPosition
    const triFace tri(currentTetIndices().faceTriIs(mesh_));
    vector n = tri.normal(mesh_.points());
    //if (mag(n) < sqr(SMALL))
    //{
    //    FatalErrorIn("wallBoundedParticle::trackFaceTri(..)")
    //        << "Small triangle." //<< static_cast<const particle&>(*this)
    //        << info()
    //        << "n:" << n
    //        << abort(FatalError);
    //}
    n /= mag(n)+VSMALL;

    // Check which edge intersects the trajectory.
    // Project trajectory onto triangle
    minEdgeI = -1;
    scalar minS = 1;        // end position

    //const point oldPosition(position());


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

        //if (mag(edgeNormal) < SMALL)
        //{
        //    FatalErrorIn("wallBoundedParticle::trackFaceTri(..)")
        //        << "Edge not perpendicular to triangle."
        //        //<< static_cast<const particle&>(*this)
        //        << info()
        //        << "triangle n:" << n
        //        << " edgeNormal:" << edgeNormal
        //        << " on tri:" << tri
        //        << " at:" << pt0
        //        << " at:" << pt1
        //        << abort(FatalError);
        //}


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


    //Pout<< "    tracked from:" << oldPosition << " to:" << position()
    //    << " projectedEnd:" << endPosition
    //    << " at s:" << minS << endl;
    //if (minEdgeI != -1)
    //{
    //    Pout<< "    on edge:" << minEdgeI
    //        << " on edge:"
    //        << mesh_.points()[tri[minEdgeI]]
    //        << mesh_.points()[tri[tri.fcIndex(minEdgeI)]]
    //        << endl;
    //}

    return minS;
}


// See if the current triangle has got a point on the
// correct side of the edge.
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
{
    //checkInside();

    //if (meshEdgeStart_ != -1 || diagEdge_ != -1)
    //{
    //    checkOnEdge();
    //}

    // Unfortunately have no access to trackdata so cannot check if particle
    // is on a wallPatch or has an mesh edge set (either of which is
    // a requirement).
}


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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallBoundedParticle::write(Ostream& os, bool writeFields) const
{
    const particle& p = static_cast<const particle&>(*this);

    if (os.format() == IOstream::ASCII)
    {
        // Write base particle
        p.write(os, writeFields);

        if (writeFields)
        {
            // Write the additional entries
            os  << token::SPACE << meshEdgeStart_
                << token::SPACE << diagEdge_;
        }
    }
    else
    {
        // Write base particle
        p.write(os, writeFields);

        // Write additional entries
        if (writeFields)
        {
            os.write
            (
                reinterpret_cast<const char*>(&meshEdgeStart_),
                sizeof(meshEdgeStart_)
              + sizeof(diagEdge_)
            );
        }
    }

    // Check state of Ostream
    os.check("wallBoundedParticle::write(Ostream& os, bool) const");
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const wallBoundedParticle& p
)
{
    // Write all data
    p.write(os, true);

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
