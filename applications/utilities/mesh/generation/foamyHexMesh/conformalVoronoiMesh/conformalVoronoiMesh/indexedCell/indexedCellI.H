/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

//#include "indexedCellChecks.H"

template<class Gt, class Cb>
Foam::tetCell CGAL::indexedCell<Gt, Cb>::unsortedVertexGlobalIndices
(
    const Foam::globalIndex& globalDelaunayVertexIndices
) const
{
    Foam::tetCell tVGI;

    for (int i = 0; i < 4; i++)
    {
        Vertex_handle v = this->vertex(i);

        // Finding the global index of each Delaunay vertex
        tVGI[i] = globalDelaunayVertexIndices.toGlobal
        (
            Foam::Pstream::myProcNo(),
            v->index()
        );
    }

    return tVGI;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Gt, class Cb>
CGAL::indexedCell<Gt, Cb>::indexedCell()
:
    Cb(),
    index_(ctUnassigned),
    filterCount_(0)
{}


template<class Gt, class Cb>
CGAL::indexedCell<Gt, Cb>::indexedCell
(
    Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3
)
:
    Cb(v0, v1, v2, v3),
    index_(ctUnassigned),
    filterCount_(0)
{}


template<class Gt, class Cb>
CGAL::indexedCell<Gt, Cb>::indexedCell
(
    Vertex_handle v0,
    Vertex_handle v1,
    Vertex_handle v2,
    Vertex_handle v3,
    Cell_handle n0,
    Cell_handle n1,
    Cell_handle n2,
    Cell_handle n3
)
:
    Cb(v0, v1, v2, v3, n0, n1, n2, n3),
    index_(ctUnassigned),
    filterCount_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Gt, class Cb>
int& CGAL::indexedCell<Gt, Cb>::cellIndex()
{
    return index_;
}


template<class Gt, class Cb>
int CGAL::indexedCell<Gt, Cb>::cellIndex() const
{
    return index_;
}


#ifdef CGAL_INEXACT

    template<class Gt, class Cb>
    const Foam::point& CGAL::indexedCell<Gt, Cb>::dual()
    {
        return reinterpret_cast<const Foam::point&>(this->circumcenter());
    }

#else

    template<class Gt, class Cb>
    const Foam::point CGAL::indexedCell<Gt, Cb>::dual()
    {
        const typename Gt::Point_3& P = this->circumcenter();

        return Foam::point
        (
            CGAL::to_double(P.x()),
            CGAL::to_double(P.y()),
            CGAL::to_double(P.z())
        );
    }

#endif


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::unassigned() const
{
    return index_ == ctUnassigned;
}


template<class Gt, class Cb>
inline int& CGAL::indexedCell<Gt, Cb>::filterCount()
{
    return filterCount_;
}


template<class Gt, class Cb>
inline int CGAL::indexedCell<Gt, Cb>::filterCount() const
{
    return filterCount_;
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::real() const
{
    return
    (
        (
            this->vertex(0)->real()
         || this->vertex(1)->real()
         || this->vertex(2)->real()
         || this->vertex(3)->real()
        )
        &&
        !(
            this->vertex(0)->farPoint()
         || this->vertex(1)->farPoint()
         || this->vertex(2)->farPoint()
         || this->vertex(3)->farPoint()
        )
    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::hasFarPoint() const
{
    return
    (
        this->vertex(0)->farPoint()
     || this->vertex(1)->farPoint()
     || this->vertex(2)->farPoint()
     || this->vertex(3)->farPoint()
    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::hasReferredPoint() const
{
    return
    (
        this->vertex(0)->referred()
     || this->vertex(1)->referred()
     || this->vertex(2)->referred()
     || this->vertex(3)->referred()
    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::hasFeaturePoint() const
{
    return
    (
        this->vertex(0)->featurePoint()
     || this->vertex(1)->featurePoint()
     || this->vertex(2)->featurePoint()
     || this->vertex(3)->featurePoint()
    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::hasSeedPoint() const
{
    return
    (
        this->vertex(0)->seedPoint()
     || this->vertex(1)->seedPoint()
     || this->vertex(2)->seedPoint()
     || this->vertex(3)->seedPoint()
    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::hasInternalPoint() const
{
    return
    (
        this->vertex(0)->internalPoint()
     || this->vertex(1)->internalPoint()
     || this->vertex(2)->internalPoint()
     || this->vertex(3)->internalPoint()
    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::hasBoundaryPoint() const
{
    return
    (
        this->vertex(0)->boundaryPoint()
     || this->vertex(1)->boundaryPoint()
     || this->vertex(2)->boundaryPoint()
     || this->vertex(3)->boundaryPoint()
    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::hasConstrainedPoint() const
{
    return
    (
        this->vertex(0)->constrained()
     || this->vertex(1)->constrained()
     || this->vertex(2)->constrained()
     || this->vertex(3)->constrained()
    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::parallelDualVertex() const
{
    return
    (
        !this->hasFarPoint()
    &&
        (
            this->vertex(0)->referred()
         || this->vertex(1)->referred()
         || this->vertex(2)->referred()
         || this->vertex(3)->referred()
        )
    &&
        (
            this->vertex(0)->real()
         || this->vertex(1)->real()
         || this->vertex(2)->real()
         || this->vertex(3)->real()
        )
    );
}


template<class Gt, class Cb>
inline Foam::label CGAL::indexedCell<Gt, Cb>::vertexLowestProc() const
{
    Foam::label lowestProc = -1;

    for (int i = 0; i < 4; ++i)
    {
        if (this->vertex(i)->referred())
        {
            lowestProc = min(lowestProc, this->vertex(i)->procIndex());
        }
    }

    return lowestProc;
}


template<class Gt, class Cb>
inline Foam::tetCell CGAL::indexedCell<Gt, Cb>::vertexGlobalIndices
(
    const Foam::globalIndex& globalDelaunayVertexIndices
) const
{
    // tetVertexGlobalIndices
    Foam::tetCell tVGI
        = unsortedVertexGlobalIndices(globalDelaunayVertexIndices);

    // bubble sort
    for (int i = 0; i < tVGI.size(); i++)
    {
        for (int j = tVGI.size() - 1 ; j > i; j--)
        {
            if (tVGI[j - 1] > tVGI[j])
            {
                Foam::Swap(tVGI[j - 1], tVGI[j]);
            }
        }
    }

    return tVGI;
}


template<class Gt, class Cb>
inline Foam::FixedList<Foam::label, 4>
CGAL::indexedCell<Gt, Cb>::globallyOrderedCellVertices
(
    const Foam::globalIndex& globalDelaunayVertexIndices
) const
{
    // tetVertexGlobalIndices
    Foam::tetCell tVGI
        = unsortedVertexGlobalIndices(globalDelaunayVertexIndices);

    Foam::FixedList<Foam::label, 4> vertexMap(Foam::identity(4));

    // bubble sort
    for (int i = 0; i < tVGI.size(); i++)
    {
        for (int j = tVGI.size() - 1 ; j > i; j--)
        {
            if (tVGI[j - 1] > tVGI[j])
            {
                Foam::Swap(tVGI[j - 1], tVGI[j]);
                Foam::Swap(vertexMap[j - 1], vertexMap[j]);
            }
        }
    }

    for (int i = 0; i < 4; i++)
    {
        tVGI[i] = vertexMap[i];
    }

    return tVGI;
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::internalOrBoundaryDualVertex() const
{
    return
    (
        this->vertex(0)->internalOrBoundaryPoint()
     || this->vertex(1)->internalOrBoundaryPoint()
     || this->vertex(2)->internalOrBoundaryPoint()
     || this->vertex(3)->internalOrBoundaryPoint()
    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::anyInternalOrBoundaryDualVertex() const
{
    return
    (
        this->vertex(0)->internalOrBoundaryPoint()
     || this->vertex(0)->externalBoundaryPoint()
     || this->vertex(1)->internalOrBoundaryPoint()
     || this->vertex(1)->externalBoundaryPoint()
     || this->vertex(2)->internalOrBoundaryPoint()
     || this->vertex(2)->externalBoundaryPoint()
     || this->vertex(3)->internalOrBoundaryPoint()
     || this->vertex(3)->externalBoundaryPoint()
    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::boundaryDualVertex() const
{
//    return
//    (
//        this->vertex(0)->boundaryPoint()
//     && this->vertex(1)->boundaryPoint()
//     && this->vertex(2)->boundaryPoint()
//     && this->vertex(3)->boundaryPoint()
//    );

    return
    (
        (
           this->vertex(0)->internalBoundaryPoint()
        || this->vertex(1)->internalBoundaryPoint()
        || this->vertex(2)->internalBoundaryPoint()
        || this->vertex(3)->internalBoundaryPoint()
        )
     && (
           this->vertex(0)->externalBoundaryPoint()
        || this->vertex(1)->externalBoundaryPoint()
        || this->vertex(2)->externalBoundaryPoint()
        || this->vertex(3)->externalBoundaryPoint()
        )
    );

//    Foam::label nBoundaryPoints = 0;
//
//    for (Foam::label i = 0; i < 4; ++i)
//    {
//        Vertex_handle v = this->vertex(i);
//
//        if (v->boundaryPoint())
//        {
//            nBoundaryPoints++;
//        }
//    }
//
//    return (nBoundaryPoints > 1);
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::baffleSurfaceDualVertex() const
{
    return
    (
        (
           this->vertex(0)->internalBaffleSurfacePoint()
        || this->vertex(1)->internalBaffleSurfacePoint()
        || this->vertex(2)->internalBaffleSurfacePoint()
        || this->vertex(3)->internalBaffleSurfacePoint()
        )
     && (
           this->vertex(0)->externalBaffleSurfacePoint()
        || this->vertex(1)->externalBaffleSurfacePoint()
        || this->vertex(2)->externalBaffleSurfacePoint()
        || this->vertex(3)->externalBaffleSurfacePoint()
        )
    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::baffleEdgeDualVertex() const
{
    return
    (
        (
           this->vertex(0)->internalBaffleEdgePoint()
        || this->vertex(1)->internalBaffleEdgePoint()
        || this->vertex(2)->internalBaffleEdgePoint()
        || this->vertex(3)->internalBaffleEdgePoint()
        )
     && (
           this->vertex(0)->externalBaffleEdgePoint()
        || this->vertex(1)->externalBaffleEdgePoint()
        || this->vertex(2)->externalBaffleEdgePoint()
        || this->vertex(3)->externalBaffleEdgePoint()
        )
    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::featureEdgeDualVertex() const
{
    return
    (
        this->vertex(0)->featureEdgePoint()
     && this->vertex(1)->featureEdgePoint()
     && this->vertex(2)->featureEdgePoint()
     && this->vertex(3)->featureEdgePoint()
    );
//        (
//            this->vertex(0)->featureEdgePoint()
//         || this->vertex(1)->featureEdgePoint()
//         || this->vertex(2)->featureEdgePoint()
//         || this->vertex(3)->featureEdgePoint()
//        )
//     && (
//            (
//                this->vertex(0)->featureEdgePoint()
//             || this->vertex(0)->featurePoint()
//            )
//         && (
//                this->vertex(1)->featureEdgePoint()
//             || this->vertex(1)->featurePoint()
//            )
//         && (
//                this->vertex(2)->featureEdgePoint()
//             || this->vertex(2)->featurePoint()
//            )
//         && (
//                this->vertex(3)->featureEdgePoint()
//             || this->vertex(3)->featurePoint()
//            )
//        )
//    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::featurePointDualVertex() const
{
    return
    (
        this->vertex(0)->featurePoint()
     && this->vertex(1)->featurePoint()
     && this->vertex(2)->featurePoint()
     && this->vertex(3)->featurePoint()
    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::nearProcBoundary() const
{
    return
    (
        this->vertex(0)->nearProcBoundary()
     || this->vertex(1)->nearProcBoundary()
     || this->vertex(2)->nearProcBoundary()
     || this->vertex(3)->nearProcBoundary()
    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::potentialCoplanarCell() const
{
    Foam::label nMasters = 0;
    Foam::label nSlaves = 0;

    Vertex_handle vM[2];
    Vertex_handle vS[2];

    for (Foam::label i = 0; i < 4; ++i)
    {
        Vertex_handle v = this->vertex(i);

        if (v->internalBoundaryPoint())
        {
            vM[nMasters] = v;
            nMasters++;
        }

        if (v->externalBoundaryPoint())
        {
            vS[nSlaves] = v;
            nSlaves++;
        }
    }

    Foam::label nPairs = 0;

    if (nMasters == 2 && nSlaves == 2)
    {
        Foam::vector vp0(Foam::vector::zero);
        Foam::vector vp1(Foam::vector::zero);

        if
        (
            vM[0]->type() == vS[0]->index()
         && vM[0]->index() == vS[0]->type()
        )
        {
            vp0 = reinterpret_cast<const Foam::point&>(vM[0]->point())
                - reinterpret_cast<const Foam::point&>(vS[0]->point());
            nPairs++;
        }
        else if
        (
            vM[0]->type() == vS[1]->index()
         && vM[0]->index() == vS[1]->type()
        )
        {
            vp0 = reinterpret_cast<const Foam::point&>(vM[0]->point())
                - reinterpret_cast<const Foam::point&>(vS[1]->point());
            nPairs++;
        }

        if
        (
            vM[1]->type() == vS[0]->index()
         && vM[1]->index() == vS[0]->type()
        )
        {
            vp1 = reinterpret_cast<const Foam::point&>(vM[1]->point())
                - reinterpret_cast<const Foam::point&>(vS[0]->point());
            nPairs++;
        }
        else if
        (
            vM[1]->type() == vS[1]->index()
         && vM[1]->index() == vS[1]->type()
        )
        {
            vp1 = reinterpret_cast<const Foam::point&>(vM[1]->point())
                - reinterpret_cast<const Foam::point&>(vS[1]->point());
            nPairs++;
        }

        if (nPairs == 2)
        {
            if (Foam::vectorTools::areParallel(vp0, vp1))
            {
                Foam::Pout<< "PARALLEL" << Foam::endl;

                return true;
            }
        }
    }

    return false;
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::featurePointExternalCell() const
{
    int featureVertex = -1;
    for (int i = 0; i < 4; ++i)
    {
        if (this->vertex(i)->constrained())
        {
            featureVertex = i;
        }
    }

    // Pick cell with a face attached to an infinite cell

    if (featureVertex != -1)
    {
        Vertex_handle v1 =
            this->vertex(Tds::vertex_triple_index(featureVertex, 0));
        Vertex_handle v2 =
            this->vertex(Tds::vertex_triple_index(featureVertex, 1));
        Vertex_handle v3 =
            this->vertex(Tds::vertex_triple_index(featureVertex, 2));

        if (v1->internalBoundaryPoint())
        {
            if
            (
                v2->externalBoundaryPoint()
             && v3->externalBoundaryPoint()
            )
            {
                return true;
            }
        }
        else if (v2->internalBoundaryPoint())
        {
            if
            (
                v1->externalBoundaryPoint()
             && v3->externalBoundaryPoint()
            )
            {
                return true;
            }
        }
        else if (v3->internalBoundaryPoint())
        {
            if
            (
                v1->externalBoundaryPoint()
             && v2->externalBoundaryPoint()
            )
            {
                return true;
            }
        }
    }

    return false;
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::featurePointInternalCell() const
{
    int featureVertex = -1;
    for (int i = 0; i < 4; ++i)
    {
        if (this->vertex(i)->constrained())
        {
            featureVertex = i;
        }
    }

    // Pick cell with a face attached to an infinite cell

    if (featureVertex != -1)
    {
        Vertex_handle v1 =
            this->vertex(Tds::vertex_triple_index(featureVertex, 0));
        Vertex_handle v2 =
            this->vertex(Tds::vertex_triple_index(featureVertex, 1));
        Vertex_handle v3 =
            this->vertex(Tds::vertex_triple_index(featureVertex, 2));

        if (v1->externalBoundaryPoint())
        {
            if
            (
                v2->internalBoundaryPoint()
             && v3->internalBoundaryPoint()
            )
            {
                return true;
            }
        }
        else if (v2->externalBoundaryPoint())
        {
            if
            (
                v1->internalBoundaryPoint()
             && v3->internalBoundaryPoint()
            )
            {
                return true;
            }
        }
        else if (v3->externalBoundaryPoint())
        {
            if
            (
                v1->internalBoundaryPoint()
             && v2->internalBoundaryPoint()
            )
            {
                return true;
            }
        }
    }

    return false;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //
