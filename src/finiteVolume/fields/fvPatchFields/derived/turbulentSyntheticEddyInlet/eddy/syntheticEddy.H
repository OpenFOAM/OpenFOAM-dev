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
    Foam::syntheticEddy

Description
    Single eddy used by turbulentSyntheticEddyInletFvPatchVectorField.

SourceFiles
    syntheticEddyI.H
    syntheticEddy.C

\*---------------------------------------------------------------------------*/

#ifndef syntheticEddy_H
#define syntheticEddy_H

#include "vector.H"
#include "point.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class syntheticEddy;

bool operator==(const syntheticEddy& a, const syntheticEddy& b);
bool operator!=(const syntheticEddy& a, const syntheticEddy& b);

/*---------------------------------------------------------------------------*\
                        Class syntheticEddy Declaration
\*---------------------------------------------------------------------------*/

class syntheticEddy
{
    // Private Data

        //- Reference position on the patch
        point position0_;

        //- Streamwise distance from the reference position
        scalar x_;

        //- Length scales in the principal frame
        vector sigma_;

        //- Intensities in the principal frame
        vector alpha_;

        //- Principal-to-global rotation tensor
        tensor Rpg_;

        //- Kernel normalisation factor
        scalar c1_;

        //- Per-eddy streamwise convection velocity. Set at spawn time. The
        //  BC's convectEddies() advances x_ by deltaT*Uconv_ each step.
        scalar Uconv_;


public:

    // Constructors

        //- Construct null
        syntheticEddy();

        //- Construct from components
        syntheticEddy
        (
            const point& position0,
            scalar x,
            const vector& sigma,
            const vector& alpha,
            const tensor& Rpg,
            scalar c1,
            scalar Uconv
        );


    // Member Functions

        //- Return the reference position
        inline const point& position0() const noexcept;

        //- Return the streamwise distance
        inline scalar x() const noexcept;

        //- Return the convection velocity
        inline scalar Uconv() const noexcept;

        //- Return the current position
        inline point position(const vector& n) const;

        //- Return the eddy volume
        inline scalar volume() const;

        //- Advance the eddy by dx
        inline void move(scalar dx);

        //- Return the velocity fluctuation at a point
        vector uPrime(const point& xp, const vector& n) const;


    // Friend Operators

        friend bool operator==(const syntheticEddy& a, const syntheticEddy& b)
        {
            return
                a.position0_ == b.position0_
             && a.x_ == b.x_
             && a.sigma_ == b.sigma_
             && a.alpha_ == b.alpha_
             && a.Rpg_ == b.Rpg_
             && a.c1_ == b.c1_
             && a.Uconv_ == b.Uconv_;
        }

        friend bool operator!=(const syntheticEddy& a, const syntheticEddy& b)
        {
            return !(a == b);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "syntheticEddyI.H"

#endif

// ************************************************************************* //
