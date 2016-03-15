/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::SpatialVector

Description
    Templated 3D spatial vector derived from VectorSpace used to represent the
    anglular and linear components of position, velocity and acceleration of
    rigid bodies.

    Reference:
    \verbatim
        Featherstone, R. (2008).
        Rigid body dynamics algorithms.
        Springer.
    \endverbatim

SourceFiles
    SpatialVectorI.H

SeeAlso
    Foam::VectorSpace
    Foam::Vector

\*---------------------------------------------------------------------------*/

#ifndef SpatialVector_H
#define SpatialVector_H

#include "Vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class SpatialVector Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class SpatialVector
:
    public VectorSpace<SpatialVector<Cmpt>, Cmpt, 6>
{

public:

    // Constructors

        //- Construct null
        inline SpatialVector();

        //- Construct initialized to zero
        inline explicit SpatialVector(const Foam::zero);

        //- Construct given VectorSpace of the same rank
        inline SpatialVector(const typename SpatialVector::vsType&);

        //- Construct from the angular and linear vector components
        inline SpatialVector
        (
            const Vector<Cmpt>& angular,
            const Vector<Cmpt>& linear
        );

        //- Construct given 6 components
        inline SpatialVector
        (
            const Cmpt& v0,
            const Cmpt& v1,
            const Cmpt& v2,
            const Cmpt& v3,
            const Cmpt& v4,
            const Cmpt& v5
        );

        //- Construct from Istream
        inline SpatialVector(Istream&);


    // Member Functions

        //- Return the angular part of the spatial vector as a vector
        inline Vector<Cmpt> angular() const;

        //- Return the linear part of the spatial vector as a vector
        inline Vector<Cmpt> linear() const;


    // Member Operators

        inline void operator=(const Foam::zero);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "SpatialVectorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
