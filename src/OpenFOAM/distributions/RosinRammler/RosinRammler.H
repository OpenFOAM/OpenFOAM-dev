/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::distributions::RosinRammler

Description
    Rosin-Rammler distribution, scaled so that it spans between a specified
    minimum and maximum value, rather than from zero to infinity

    \f[
        PDF(x) = \left( \frac{n}{d} \right) \left( \frac{x}{d} \right)^{n-1} \\
        \exp \left( - \left(\frac{x}{d} \right)^n \right)
    \f]
    \f[
        CDF(x) = 1 - \exp \left( - \left( \frac{x}{d} \right)^n \right)
    \f]

Usage
    Example usage:
    \verbatim
    {
        type    RosinRammler;
        Q       0;
        min     0.00001;
        max     0.00015;
        d       0.00014;
        n       2;
    }
    \endverbatim

SourceFiles
    RosinRammler.C

See also
    Foam::distribution

\*---------------------------------------------------------------------------*/

#ifndef RosinRammler_H
#define RosinRammler_H

#include "unintegrable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{

/*---------------------------------------------------------------------------*\
                        Class RosinRammler Declaration
\*---------------------------------------------------------------------------*/

class RosinRammler
:
    public FieldDistribution<unintegrableForNonZeroQ, RosinRammler>
{
    // Private Data

        //- Minimum value
        const scalar min_;

        //- Maximum value
        const scalar max_;

        //- Scale parameter
        const scalar d_;

        //- Shape parameter
        const scalar n_;


protected:

    // Protected Member Functions

        //- Return values of the un-normalised PDF for the given size exponent
        //  and x-coordinates.
        virtual tmp<scalarField> phi
        (
            const label q,
            const scalarField& x
        ) const;

        //- Return values of the un-normalised CDF for the given size exponent
        //  and x-coordinates.
        virtual tmp<scalarField> Phi
        (
            const label q,
            const scalarField& x
        ) const;


public:

    //- Runtime type information
    TypeName("RosinRammler");


    // Constructors

        //- Construct from a dictionary
        RosinRammler
        (
            const dictionary& dict,
            randomGenerator& rndGen,
            const label sampleQ
        );

        //- Construct copy
        RosinRammler(const RosinRammler& d, const label sampleQ);

        //- Construct and return a clone
        virtual autoPtr<distribution> clone(const label sampleQ) const
        {
            return autoPtr<distribution>(new RosinRammler(*this, sampleQ));
        }


    //- Destructor
    virtual ~RosinRammler();


    // Member Functions

        //- Sample the distribution
        virtual scalar sample() const;

        //- Sample the distribution
        using FieldDistribution<unintegrableForNonZeroQ, RosinRammler>::sample;

        //- Return the minimum value
        virtual scalar min() const;

        //- Return the maximum value
        virtual scalar max() const;

        //- Return coordinates to plot across the range of the distribution
        virtual tmp<scalarField> x(const label n) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace distributions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
