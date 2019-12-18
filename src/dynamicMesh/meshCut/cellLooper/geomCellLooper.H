/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::geomCellLooper

Description
    Implementation of cellLooper. Does pure geometric cut through cell.

    Handles all cell shapes in the same way: cut edges with plane through
    cell centre and normal in direction of provided direction. Snaps cuts
    close to edge endpoints (close = snapTol * minEdgeLen) to vertices.

    Currently determines cuts through edges (and edgeendpoints close to plane)
    in random order and then sorts them acc. to angle. Could be converted to
    use walk but problem is that face can be cut multiple times (since does
    not need to be convex). Another problem is that edges parallel to plane
    might not be cut. So these are handled by looking at the distance from
    edge endpoints to the plane.

SourceFiles
    geomCellLooper.C

\*---------------------------------------------------------------------------*/

#ifndef geomCellLooper_H
#define geomCellLooper_H

#include "cellLooper.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class plane;

/*---------------------------------------------------------------------------*\
                           Class geomCellLooper Declaration
\*---------------------------------------------------------------------------*/

class geomCellLooper
:
    public cellLooper
{
    // Static

        //- Tolerance for point equal test. Fraction of edge length.
        static const scalar pointEqualTol_;

        //- Tolerance for cut through edges to get snapped to edge end point.
        //  Fraction of length of minimum connected edge length.
        static scalar snapTol_;


    // Private Member Functions

        //- Min length of attached edges
        scalar minEdgeLen(const label vertI) const;

        //- Return true and set weight if edge is cut
        bool cutEdge
        (
            const plane& cutPlane,
            const label edgeI,
            scalar& weight
        ) const;

        //- Snaps cut through edge by cut through vertex (if weight closer than
        //  tol to 0 or 1). Returns vertex label snapped to or -1.
        label snapToVert
        (
            const scalar tol,
            const label edgeI,
            const scalar weight
        ) const;

        //- Gets two (random) vectors perpendicular to n and each other to be
        //  used as base.
        void getBase
        (
            const vector& n,
            vector& e0,
            vector& e1
        ) const;

        //- Return true if the cut edge at loop[index] is in between the cuts
        //  through the edge end points.
        bool edgeEndsCut(const labelList&, const label index) const;


public:

    //- Runtime type information
    TypeName("geomCellLooper");


    // Static Functions

        static scalar snapTol()
        {
            return snapTol_;
        }

        static void setSnapTol(const scalar tol)
        {
            snapTol_ = tol;
        }


    // Constructors

        //- Construct from components
        geomCellLooper(const polyMesh& mesh);

        //- Disallow default bitwise copy construction
        geomCellLooper(const geomCellLooper&) = delete;


    //- Destructor
    virtual ~geomCellLooper();


    // Member Functions

        //- Create cut along circumference of celli. Gets current mesh cuts.
        //  Cut along circumference is expressed as loop of cuts plus weights
        //  for cuts along edges (only valid for edge cuts).
        //  Return true if successful cut.
        virtual bool cut
        (
            const vector& refDir,
            const label celli,
            const boolList& vertIsCut,
            const boolList& edgeIsCut,
            const scalarField& edgeWeight,

            labelList& loop,
            scalarField& loopWeights
        ) const;

        //- Same but now also base point of cut provided (instead of always
        //  cell centre)
        virtual bool cut
        (
            const plane& cutPlane,
            const label celli,
            const boolList& vertIsCut,
            const boolList& edgeIsCut,
            const scalarField& edgeWeight,

            labelList& loop,
            scalarField& loopWeights
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const geomCellLooper&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
