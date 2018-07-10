/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

Class
    CGAL::indexedCell

Description
    An indexed form of CGAL::Triangulation_cell_base_3<K> used to keep
    track of the Delaunay cells (tets) in the tessellation.

SourceFiles
    indexedCellI.H
    indexedCell.C

\*---------------------------------------------------------------------------*/

#ifndef indexedCell_H
#define indexedCell_H

#include "CGAL/Triangulation_3.h"
#include "CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h"
#include "indexedVertex.H"
#include "List.H"
#include "globalIndex.H"
#include "Pstream.H"
#include "Swap.H"
#include "InfoProxy.H"
#include "tetCell.H"
#include "typeInfo.H"
#include "vectorTools.H"
#include "indexedCellEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace CGAL
{
template<class Gt, class Cb> class indexedCell;
}

namespace Foam
{
class Ostream;
template<class Gt, class Cb> Ostream& operator<<
(
    Ostream&,
    const Foam::InfoProxy<CGAL::indexedCell<Gt, Cb>>&
);
}


namespace CGAL
{

/*---------------------------------------------------------------------------*\
                          Class indexedCell Declaration
\*---------------------------------------------------------------------------*/

template
<
    class Gt,
    class Cb=CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3<Gt>
>
class indexedCell
:
    public Foam::indexedCellEnum,
    public Cb
{
    // Private data

        //- The index for this Delaunay tetrahedral cell.  Type information is
        //- Also carried:

        //  ctFar : the dual point of this cell does not form part of the
        //          internal or boundary of the dual mesh
        //  >=0   : the (local) index of an internal or boundary dual point,
        //           not on a processor face
        //  < 0 && > ctFar : the (global) index of a dual point on a processor
        //                   face

        Foam::label index_;

        //- The number of times that this Delaunay cell has been limited
        //  during filtering
        int filterCount_;


    // Private member functions

        //- Same as globallyOrderedCellVertices but without sorting
        Foam::tetCell unsortedVertexGlobalIndices
        (
            const Foam::globalIndex& globalDelaunayVertexIndices
        ) const;


public:

    typedef typename Cb::Triangulation_data_structure    Tds;
    typedef typename Cb::Vertex_handle                   Vertex_handle;
    typedef typename Cb::Cell_handle                     Cell_handle;

    template<class TDS2>
    struct Rebind_TDS
    {
        typedef typename Cb::template Rebind_TDS<TDS2>::Other Cb2;
        typedef indexedCell<Gt, Cb2>  Other;
    };


    // Constructors

        inline indexedCell();

        inline indexedCell
        (
            Vertex_handle v0,
            Vertex_handle v1,
            Vertex_handle v2,
            Vertex_handle v3
        );

        inline indexedCell
        (
            Vertex_handle v0,
            Vertex_handle v1,
            Vertex_handle v2,
            Vertex_handle v3,
            Cell_handle n0,
            Cell_handle n1,
            Cell_handle n2,
            Cell_handle n3
        );


    // Member Functions

        inline Foam::label& cellIndex();

        inline Foam::label cellIndex() const;

#ifdef CGAL_INEXACT
        inline const Foam::point& dual();
#else
        inline const Foam::point dual();
#endif

        inline bool unassigned() const;

        inline int& filterCount();

        inline int filterCount() const;

        //- Is the Delaunay cell real, i.e. any real vertex
        inline bool real() const;

        //- Does the Delaunay cell have a far point
        inline bool hasFarPoint() const;

        //- Does the Delaunay cell have a referred point
        inline bool hasReferredPoint() const;

        //- Does the Delaunay cell have a feature point
        inline bool hasFeaturePoint() const;

        //- Does the Delaunay cell have a seed point
        inline bool hasSeedPoint() const;

        inline bool hasInternalPoint() const;
        inline bool hasBoundaryPoint() const;

        inline bool hasConstrainedPoint() const;

        //- Does the Dual vertex form part of a processor patch
        inline bool parallelDualVertex() const;

        inline Foam::label vertexLowestProc() const;

        //- Using the globalIndex object, return a list of four (sorted) global
        //  Delaunay vertex indices that uniquely identify this tet in parallel
        inline Foam::tetCell vertexGlobalIndices
        (
            const Foam::globalIndex& globalDelaunayVertexIndices
        ) const;

        //- Using the globalIndex object, return a list of four vertices with
        //  so that the cell has a consistent orientation in parallel.
        inline Foam::FixedList<Foam::label, 4> globallyOrderedCellVertices
        (
            const Foam::globalIndex& globalDelaunayVertexIndices
        ) const;

        //- Is the Delaunay cell part of the final dual mesh, i.e. any vertex
        //  form part of the internal or boundary definition
        inline bool internalOrBoundaryDualVertex() const;

        //- Is the Delaunay cell real or referred (or mixed), i.e. all vertices
        //  form part of the real or referred internal or boundary definition
        inline bool anyInternalOrBoundaryDualVertex() const;

        //- A dual vertex on the boundary will result from a Delaunay cell with
        //  least one Delaunay vertex outside and at least one inside
        inline bool boundaryDualVertex() const;

        inline bool baffleSurfaceDualVertex() const;

        inline bool baffleEdgeDualVertex() const;

        //- A dual vertex on a feature edge will result from this Delaunay cell
        inline bool featureEdgeDualVertex() const;

        //- A dual vertex on a feature point will result from this Delaunay cell
        inline bool featurePointDualVertex() const;

        inline bool nearProcBoundary() const;

        inline bool potentialCoplanarCell() const;

        inline bool featurePointExternalCell() const;
        inline bool featurePointInternalCell() const;


    // Info

        //- Return info proxy.
        //  Used to print indexedCell information to a stream
        Foam::InfoProxy<indexedCell<Gt, Cb>> info() const
        {
            return *this;
        }

        friend Foam::Ostream& Foam::operator<< <Gt, Cb>
        (
            Foam::Ostream&,
            const Foam::InfoProxy<indexedCell<Gt, Cb>>&
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CGAL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "indexedCellI.H"

#ifdef NoRepository
    #include "indexedCell.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
