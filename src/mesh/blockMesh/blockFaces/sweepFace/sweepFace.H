/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::blockFaces::sweepFace

Description
    This curved face sweeps the points of two of the face's four edges over the
    surface of the face, using the other two edges to define the
    transformation. It does a better job than the standard interpolation
    algorithm of conforming to cylinders, cones, toroids and similar.

    The algorithm is asymmetric, and depends on the order in which the vertices
    are supplied. E.g., if vertices (0 1 2 3) are supplied, then that means
    edge (0 1) will be swept towards edge (2 3) in the direction defined by
    edges (0 3) and (1 2).

    The sweeping is controlled by a single scalar parameter, the rotation
    fraction, which must take a value between 0 and 1. This fraction determines
    to what extent the swept edges are rotated around their own axes.

SourceFiles
    sweepFace.C

\*---------------------------------------------------------------------------*/

#ifndef sweepFace_blockFace_H
#define sweepFace_blockFace_H

#include "blockFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockFaces
{

/*---------------------------------------------------------------------------*\
                         Class sweepFace Declaration
\*---------------------------------------------------------------------------*/

class sweepFace
:
    public blockFace
{
    // Private Data

        //- Rotation fraction
        const scalar rotationFraction_;


public:

    //- Runtime type information
    TypeName("sweep");


    // Constructors

        //- Construct from Istream setting pointsList
        sweepFace
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaceList& geometry,
            Istream&
        );

        //- Disallow default bitwise copy construction
        sweepFace(const sweepFace&) = delete;


    //- Destructor
    virtual ~sweepFace()
    {}


    // Member Functions

        //- Project the given points onto the surface
        virtual void project
        (
            const blockDescriptor&,
            const label blockFacei,
            pointField& points
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sweepFace&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockFaces
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
