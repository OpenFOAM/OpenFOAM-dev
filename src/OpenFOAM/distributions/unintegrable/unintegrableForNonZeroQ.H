/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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
    Foam::distributions::unintegrableForNonZeroQ

Description
    Base class for distributions that have a closed integral form for the
    cumulative density function (CDF) when the effective size exponent is zero,
    but not otherwise.

SourceFiles
    unintegrableForNonZeroQ.C

See also
    Foam::distribution

\*---------------------------------------------------------------------------*/

#ifndef unintegrableForNonZeroQ_H
#define unintegrableForNonZeroQ_H

#include "unintegrable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{

/*---------------------------------------------------------------------------*\
                        Class unintegrableForNonZeroQ Declaration
\*---------------------------------------------------------------------------*/

class unintegrableForNonZeroQ
:
    public unintegrable
{
protected:

    // Protected Member Functions

        //- Return values of the un-normalised CDF for the given size exponent
        //  and x-coordinates.
        virtual tmp<scalarField> Phi
        (
            const label q,
            const scalarField& x
        ) const;

        //- Return values of the un-normalised CDF for zero effective size
        //  exponent and given x-coordinates
        virtual tmp<scalarField> PhiForZeroQ(const scalarField& x) const = 0;

        //- Return values of the un-normalised CDF at the minimum and maximum
        //  x-coordinates for the given size exponent.
        virtual Pair<scalar> Phi01(const label q) const;

        //- Access cached values of the un-normalised CDF at the minimum and
        //  maximum x-coordinates.
        using unintegrable::Phi01;


public:

    // Constructors

        //- Inherit constructors
        using unintegrable::unintegrable;


    // Member Functions

        //- Sample the distribution
        virtual scalar sample() const;

        //- Sample the distribution
        virtual scalar sampleForZeroQ() const = 0;

        //- Return the cumulative density function at the given coordinates
        virtual tmp<scalarField> CDF(const scalarField& x) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace distributions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
