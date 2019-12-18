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
    Foam::blockEdges::arcEdge

Description
    An arcEdge between two points on a circle. The arc is defined either by a
    third point that the arc passes through, or by the angle of the sector and
    the axis of the circle.

    If the angle-axis specification is used, the axis may not be not
    perpendicular to the line connecting the end points. In that case, the axis
    is considered to be that of a cylinder, and the arc represents a portion of
    a helix on the surface of that cylinder.

SourceFiles
    arcEdge.C

\*---------------------------------------------------------------------------*/

#ifndef blockEdges_arcEdge_H
#define blockEdges_arcEdge_H

#include "blockEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{

/*---------------------------------------------------------------------------*\
                           Class arcEdge Declaration
\*---------------------------------------------------------------------------*/

class arcEdge
:
    public blockEdge
{
    // Private Data

        //- The start point of the arc
        const point p0_;

        //- The end point of the arc
        const point p1_;

        //- The centre of the circle of which the arc is a part
        point centre_;

        //- The axis of the circle of which the arc is a part
        vector axis_;

        //- The sector angle of the arc
        scalar theta_;

        //- The axial length of the arc/helix
        scalar length_;


    // Private Member Functions

        //- Calculate the coordinate system, angle and radius from a mid point
        //  somewhere within the arc
        void calc(const point& pM);

        //- Calculate the coordinate system, angle and radius from a sector
        //  angle and a normal to the plane of the arc
        void calc(const scalar theta, const vector& axis);


public:

    //- Runtime type information
    TypeName("arc");


    // Constructors

        //- Construct from Istream setting pointsList
        arcEdge
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces& geometry,
            const pointField& points,
            Istream&
        );

        //- Disallow default bitwise copy construction
        arcEdge(const arcEdge&) = delete;


    //- Destructor
    virtual ~arcEdge()
    {}


    // Member Functions

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        point position(const scalar) const;

        //- Return the length of the curve
        scalar length() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const arcEdge&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End of namespace blockEdges
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
