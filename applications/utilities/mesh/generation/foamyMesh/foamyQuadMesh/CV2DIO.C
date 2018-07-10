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
#include "OFstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CV2D::writePoints(const fileName& fName, bool internalOnly) const
{
    Info<< "Writing points to " << fName << nl << endl;
    OFstream str(fName);

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (!internalOnly || vit->internalOrBoundaryPoint())
        {
            meshTools::writeOBJ(str, toPoint3D(vit->point()));
        }
    }
}


void Foam::CV2D::writeTriangles(const fileName& fName, bool internalOnly) const
{
    Info<< "Writing triangles to " << fName << nl << endl;
    OFstream str(fName);

    labelList vertexMap(number_of_vertices(), -2);
    label verti = 0;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (!internalOnly || !vit->farPoint())
        {
            vertexMap[vit->index()] = verti++;
            meshTools::writeOBJ(str, toPoint3D(vit->point()));
        }
    }

    for
    (
        Triangulation::Finite_faces_iterator fit = finite_faces_begin();
        fit != finite_faces_end();
        ++fit
    )
    {
        if
        (
            !internalOnly
         || (
                fit->vertex(0)->internalOrBoundaryPoint()
             || fit->vertex(1)->internalOrBoundaryPoint()
             || fit->vertex(2)->internalOrBoundaryPoint()
            )
        )
        {
            str << "f";
            for (label i = 0; i < 3; ++i)
            {
                str << " " << vertexMap[fit->vertex(i)->index()] + 1;
            }
            str << nl;
        }
    }
}


void Foam::CV2D::writeFaces(const fileName& fName, bool internalOnly) const
{
    Info<< "Writing dual faces to " << fName << nl << endl;
    OFstream str(fName);

    label dualVerti = 0;

    for
    (
        Triangulation::Finite_faces_iterator fit = finite_faces_begin();
        fit != finite_faces_end();
        ++fit
    )
    {
        if
        (
            !internalOnly
         || (
                fit->vertex(0)->internalOrBoundaryPoint()
             || fit->vertex(1)->internalOrBoundaryPoint()
             || fit->vertex(2)->internalOrBoundaryPoint()
            )
        )
        {
            fit->faceIndex() = dualVerti++;
            meshTools::writeOBJ(str, toPoint3D(circumcenter(fit)));
        }
        else
        {
            fit->faceIndex() = -1;
        }
    }

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (!internalOnly || vit->internalOrBoundaryPoint())
        {
            Face_circulator fcStart = incident_faces(vit);
            Face_circulator fc = fcStart;

            str<< 'f';

            do
            {
                if (!is_infinite(fc))
                {
                    if (fc->faceIndex() < 0)
                    {
                        FatalErrorInFunction
                         << "Dual face uses vertex defined by a triangle"
                            " defined by an external point"
                            << exit(FatalError);
                    }

                    str<< ' ' << fc->faceIndex() + 1;
                }
            } while (++fc != fcStart);

            str<< nl;
        }
    }
}


void Foam::CV2D::extractPatches
(
    wordList& patchNames,
    labelList& patchSizes,
    EdgeMap<label>& mapEdgesRegion,
    EdgeMap<label>& indirectPatchEdge
) const
{
    label nPatches = qSurf_.patchNames().size() + 1;
    label defaultPatchIndex = qSurf_.patchNames().size();

    patchNames.setSize(nPatches);
    patchSizes.setSize(nPatches, 0);
    mapEdgesRegion.clear();

    const wordList& existingPatches = qSurf_.patchNames();

    forAll(existingPatches, sP)
    {
        patchNames[sP] = existingPatches[sP];
    }

    patchNames[defaultPatchIndex] = "CV2D_default_patch";

    for
    (
        Triangulation::Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        Face_handle fOwner = eit->first;
        Face_handle fNeighbor = fOwner->neighbor(eit->second);

        Vertex_handle vA = fOwner->vertex(cw(eit->second));
        Vertex_handle vB = fOwner->vertex(ccw(eit->second));

        if
        (
            (vA->internalOrBoundaryPoint() && !vB->internalOrBoundaryPoint())
         || (vB->internalOrBoundaryPoint() && !vA->internalOrBoundaryPoint())
        )
        {
            point ptA = toPoint3D(vA->point());
            point ptB = toPoint3D(vB->point());

            label patchIndex = qSurf_.findPatch(ptA, ptB);

            if (patchIndex == -1)
            {
                patchIndex = defaultPatchIndex;

                WarningInFunction
                    << "Dual face found that is not on a surface "
                    << "patch. Adding to CV2D_default_patch."
                    << endl;
            }

            edge e(fOwner->faceIndex(), fNeighbor->faceIndex());
            patchSizes[patchIndex]++;
            mapEdgesRegion.insert(e, patchIndex);

            if (!pointPair(*vA, *vB))
            {
                indirectPatchEdge.insert(e, 1);
            }
        }
    }
}


void Foam::CV2D::calcDual
(
    point2DField& dualPoints,
    faceList& dualFaces,
    wordList& patchNames,
    labelList& patchSizes,
    EdgeMap<label>& mapEdgesRegion,
    EdgeMap<label>& indirectPatchEdge
) const
{
    // Dual points stored in triangle order.
    dualPoints.setSize(number_of_faces());
    label dualVerti = 0;

    for
    (
        Triangulation::Finite_faces_iterator fit = finite_faces_begin();
        fit != finite_faces_end();
        ++fit
    )
    {
        if
        (
            fit->vertex(0)->internalOrBoundaryPoint()
         || fit->vertex(1)->internalOrBoundaryPoint()
         || fit->vertex(2)->internalOrBoundaryPoint()
        )
        {
            fit->faceIndex() = dualVerti;

            dualPoints[dualVerti++] = toPoint2D(circumcenter(fit));
        }
        else
        {
            fit->faceIndex() = -1;
        }
    }

    dualPoints.setSize(dualVerti);

    extractPatches(patchNames, patchSizes, mapEdgesRegion, indirectPatchEdge);

    forAll(patchNames, patchi)
    {
        Info<< "Patch " << patchNames[patchi]
            << " has size " << patchSizes[patchi] << endl;
    }

    // Create dual faces
    // ~~~~~~~~~~~~~~~~~

    dualFaces.setSize(number_of_vertices());
    label dualFacei = 0;
    labelList faceVerts(maxNvert);

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalOrBoundaryPoint())
        {
            Face_circulator fcStart = incident_faces(vit);
            Face_circulator fc = fcStart;
            label verti = 0;

            do
            {
                if (!is_infinite(fc))
                {
                    if (fc->faceIndex() < 0)
                    {
                        FatalErrorInFunction
                         << "Dual face uses vertex defined by a triangle"
                            " defined by an external point"
                            << exit(FatalError);
                    }

                    // Look up the index of the triangle
                    faceVerts[verti++] = fc->faceIndex();
                }
            } while (++fc != fcStart);

            if (faceVerts.size() > 2)
            {
                dualFaces[dualFacei++] =
                    face(labelList::subList(faceVerts, verti));
            }
            else
            {
                Info<< "From triangle point:" << vit->index()
                    << " coord:" << toPoint2D(vit->point())
                    << " generated illegal dualFace:" << faceVerts
                    << endl;
            }
        }
    }

    dualFaces.setSize(dualFacei);
}


void Foam::CV2D::writePatch(const fileName& fName) const
{
    point2DField dual2DPoints;
    faceList dualFaces;
    wordList patchNames;
    labelList patchSizes;
    EdgeMap<label> mapEdgesRegion;
    EdgeMap<label> indirectPatchEdge;

    calcDual
    (
        dual2DPoints,
        dualFaces,
        patchNames,
        patchSizes,
        mapEdgesRegion,
        indirectPatchEdge
    );

    pointField dualPoints(dual2DPoints.size());
    forAll(dualPoints, ip)
    {
        dualPoints[ip] = toPoint3D(dual2DPoints[ip]);
    }

    // Dump as primitive patch to be read by extrudeMesh.
    OFstream str(fName);

    Info<< "Writing patch to be used with extrudeMesh to file " << fName
        << endl;

    str << dualPoints << nl << dualFaces << nl;
}


// ************************************************************************* //
