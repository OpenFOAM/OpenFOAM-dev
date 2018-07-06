/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "CV2D.H"
#include "plane.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::CV2D::on2DLine(const point2D& p, const linePointRef& line)
{
    const point2D& a = toPoint2D(line.start());
    const point2D& b = toPoint2D(line.end());

    if
    (
        p.x() < min(a.x(), b.x())
     || p.x() > max(a.x(), b.x())
     || p.y() < min(a.y(), b.y())
     || p.y() > max(a.y(), b.y())
    )
    {
        return false;
    }

    return true;
}


void Foam::CV2D::insertFeaturePoints()
{
    featurePoints_.clear();
    label nVert = number_of_vertices();

    const PtrList<extendedFeatureEdgeMesh>& feMeshes
    (
        qSurf_.features()
    );

    if (feMeshes.empty())
    {
        WarningInFunction
            << "Extended Feature Edge Mesh is empty so no feature points will "
            << "be found." << nl
            << "    Use: featureMethod extendedFeatureEdgeMesh;" << nl
            << endl;
    }

    forAll(feMeshes, i)
    {
        const extendedFeatureEdgeMesh& feMesh(feMeshes[i]);
        const edgeList& edges = feMesh.edges();
        const pointField& points = feMesh.points();

        if (debug)
        {
            label nConvex = feMesh.concaveStart() - feMesh.convexStart();
            label nConcave = feMesh.mixedStart() - feMesh.concaveStart();
            label nMixed = feMesh.nonFeatureStart() - feMesh.mixedStart();
            label nExternal = feMesh.internalStart() - feMesh.externalStart();
            label nInternal = feMesh.flatStart() - feMesh.internalStart();
            label nFlat = feMesh.openStart() - feMesh.flatStart();
            label nOpen = feMesh.multipleStart() - feMesh.openStart();
            label nMultiple = edges.size() - feMesh.multipleStart();

            Info<< "Inserting Feature Points:" << nl
                << "    Convex points: " << nConvex << nl
                << "   Concave points: " << nConcave << nl
                << "     Mixed points: " << nMixed << nl
                << "   External edges: " << nExternal << nl
                << "   Internal edges: " << nInternal << nl
                << "       Flat edges: " << nFlat << nl
                << "       Open edges: " << nOpen << nl
                << "   Multiple edges: " << nMultiple << endl;
        }

        // Args: (base point, normal)
        // TODO: allow user to input this
        plane zPlane(vector(0, 0, z_), vector(0, 0, 1));

        if (debug)
        {
            Info<< "    plane: " << zPlane << " " << z_ << endl;
        }

        forAll(edges, edgeI)
        {
            const edge& e = feMesh.edges()[edgeI];

            const point& ep0 = points[e.start()];
            const point& ep1 = points[e.end()];

            const linePointRef line(ep0, ep1);

            scalar intersect = zPlane.lineIntersect(line);

            point2D featPoint = toPoint2D(intersect * (ep1 - ep0) + ep0);

            if (on2DLine(featPoint, line))
            {
                vector2DField fpn = toPoint2D(feMesh.edgeNormals(edgeI));

                vector2D cornerNormal = sum(fpn);
                cornerNormal /= mag(cornerNormal);

                if (debug)
                {
                    Info<< nl << "    line: " << line << nl
                        << "        vec: " << line.vec() << nl
                        << "    featurePoint: " << featPoint << nl
                        << "    line length: " << line.mag() << nl
                        << "    intersect: " << intersect << endl;
                }

                if
                (
                    feMesh.getEdgeStatus(edgeI)
                 == extendedFeatureEdgeMesh::EXTERNAL
                )
                {
                    // Convex Point
                    Foam::point2D internalPt =
                        featPoint - meshControls().ppDist()*cornerNormal;

                    if (debug)
                    {
                        Info<< "PREC: " << internalPt << nl
                            << "    : " << featPoint << nl
                            << "    : " << meshControls().ppDist() << nl
                            << "    : " << cornerNormal << endl;
                    }

                    featurePoints_.push_back
                    (
                        Vb
                        (
                            toPoint(internalPt),
                            nVert,
                            nVert + 1
                        )
                    );
                    label masterPtIndex = nVert++;

                    forAll(fpn, nI)
                    {
                        const vector n3D(fpn[nI][0], fpn[nI][1], 0.0);

                        plane planeN = plane(toPoint3D(featPoint), n3D);

                        Foam::point2D externalPt =
                            internalPt
                          + (
                                2.0
                              * planeN.distance(toPoint3D(internalPt))
                              * fpn[nI]
                            );

                        featurePoints_.push_back
                        (
                            Vb
                            (
                                toPoint(externalPt),
                                nVert++,
                                masterPtIndex
                            )
                        );

                        if (debug)
                        {
                            Info<< "  side point: " << externalPt << endl;
                        }
                    }

                    if (debug)
                    {
                        Info<< "Convex Point: " << featPoint << nl
                            << "  corner norm: " << cornerNormal << nl
                            << "    reference: " << internalPt << endl;
                    }
                }
                else if
                (
                    feMesh.getEdgeStatus(edgeI)
                 == extendedFeatureEdgeMesh::INTERNAL
                )
                {
                    // Concave Point
                    Foam::point2D externalPt =
                        featPoint + meshControls().ppDist()*cornerNormal;

                    Foam::point2D refPt =
                        featPoint - meshControls().ppDist()*cornerNormal;

                    label slavePointIndex = 0;

                    scalar totalAngle =
                        radToDeg
                        (
                            constant::mathematical::pi
                          + acos(mag(fpn[0] & fpn[1]))
                        );

                    // Number of quadrants the angle should be split into
                    int nQuads =
                        int(totalAngle/meshControls().maxQuadAngle()) + 1;

                    // The number of additional master points needed to
                    // obtain the required number of quadrants.
                    int nAddPoints = min(max(nQuads - 2, 0), 2);

                    // index of reflMaster
                    label reflectedMaster = nVert + 2 + nAddPoints;

                    if (debug)
                    {
                        Info<< "Concave Point: " <<  featPoint << nl
                            << "  corner norm: " << cornerNormal << nl
                            << "     external: " << externalPt << nl
                            << "    reference: " << refPt << nl
                            << "        angle: " << totalAngle << nl
                            << "       nQuads: " << nQuads << nl
                            << "   nAddPoints: " << nAddPoints << endl;
                    }

                    forAll(fpn, nI)
                    {
                        const vector n3D(fpn[nI][0], fpn[nI][1], 0.0);

                        plane planeN = plane(toPoint3D(featPoint), n3D);

                        Foam::point2D internalPt =
                            externalPt
                          - (
                                2.0
                              * planeN.distance(toPoint3D(externalPt))
                              * fpn[nI]
                            );

                        featurePoints_.push_back
                        (
                            Vb
                            (
                                toPoint(internalPt),
                                nVert,
                                reflectedMaster
                            )
                        );
                        slavePointIndex = nVert++;

                        if (debug)
                        {
                            Info<< "Internal Point: " <<  internalPt << endl;
                        }
                    }

                    if (nAddPoints == 1)
                    {
                        // One additional point is the reflection of the slave
                        // point, i.e., the original reference point
                        featurePoints_.push_back
                        (
                            Vb
                            (
                                toPoint(refPt),
                                nVert++,
                                reflectedMaster
                            )
                        );

                        if (debug)
                        {
                           Info<< "ref Point: " <<  refPt << endl;
                        }
                    }
                    else if (nAddPoints == 2)
                    {
                       point2D reflectedAa =
                           refPt - ((featPoint - externalPt) & fpn[1])*fpn[1];

                       featurePoints_.push_back
                       (
                           Vb
                           (
                               toPoint(reflectedAa),
                               nVert++,
                               reflectedMaster
                           )
                       );

                       point2D reflectedBb =
                           refPt - ((featPoint - externalPt) & fpn[0])*fpn[0];

                       featurePoints_.push_back
                       (
                           Vb
                           (
                               toPoint(reflectedBb),
                               nVert++,
                               reflectedMaster
                           )
                       );

                       if (debug)
                       {
                           Info<< "refA Point: " <<  reflectedAa << nl
                               << "refb Point: " <<  reflectedBb << endl;
                       }
                    }

                    featurePoints_.push_back
                    (
                        Vb
                        (
                            toPoint(externalPt),
                            nVert++,
                            slavePointIndex
                        )
                    );
                }
                else
                {
                    WarningInFunction
                        << "Feature Edge " << edges[edgeI] << nl
                        << "    points(" << points[edges[edgeI].start()]
                        << ", " << points[edges[edgeI].end()] << ")" << nl
                        << "    is not labelled as either concave or convex, it"
                        << " is labelled as (#2 = flat): "
                        << feMesh.getEdgeStatus(edgeI) << endl;
                }
            }
            else
            {
                WarningInFunction
                      << "Point " << featPoint << " is not on the line "
                      << line << endl;
            }
        }
    }

    // Insert the feature points.
    reinsertFeaturePoints();

    if (meshControls().objOutput())
    {
        writePoints("feat_allPoints.obj", false);
        writeFaces("feat_allFaces.obj", false);
        writeFaces("feat_faces.obj", true);
        writeTriangles("feat_triangles.obj", true);
    }
}


void Foam::CV2D::reinsertFeaturePoints()
{
    for
    (
        std::list<Vb>::iterator vit=featurePoints_.begin();
        vit != featurePoints_.end();
        ++vit
    )
    {
        insertPoint
        (
            toPoint2D(vit->point()),
            vit->index(),
            vit->type()
        );
    }
}


// ************************************************************************* //
