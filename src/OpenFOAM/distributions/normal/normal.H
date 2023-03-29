/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::distributions::normal

Description
    Normal distribution, scaled so that it spans between a specified
    minimum and maximum value, rather than from zero to infinity

    \f[
        PDF(x) = \frac{1}{\sigma \sqrt{2 \pi}} \\
        \exp \left( \frac{1}{2} \left( \frac{x - \mu}{\sigma} \right)^2 \right)
    \f]

Usage
    Example usage:
    \verbatim
    {
        type    normal;
        Q       0;
        min     0.001;
        max     0.019;
        mu      0.011;
        sigma   0.003;
    }
    \endverbatim

SourceFiles
    normal.C

See also
    Foam::distribution

\*---------------------------------------------------------------------------*/

#ifndef normal_H
#define normal_H

#include "unintegrable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{

/*---------------------------------------------------------------------------*\
                        Class normal Declaration
\*---------------------------------------------------------------------------*/

class normal
:
    public FieldDistribution<unintegrableForNonZeroQ, normal>
{
    // Private Data

        //- Minimum value
        const scalar min_;

        //- Maximum value
        const scalar max_;

        //- Mean
        const scalar mu_;

        //- Standard deviation
        const scalar sigma_;

        //- Constant for approximate error function and inverse error function
        static const scalar a_;


    // Private Member Functions

        //- Approximate error function
        static tmp<scalarField> approxErf(const scalarField& x);

        //- Approximate error function inverse
        static scalar approxErfInv(const scalar y);

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

        //- Sample the distribution for zero effective size exponent
        scalar sampleForZeroQ(const scalar s) const;


public:

    //- Runtime type information
    TypeName("normal");


    //- Permit the multiNormal distribution to use private parts of this class
    friend class multiNormal;


    // Constructors

        //- Construct from a dictionary
        normal
        (
            const dictionary& dict,
            Random& rndGen,
            const label sampleQ
        );

        //- Construct from components
        normal
        (
            Random& rndGen,
            const label Q,
            const label sampleQ,
            const label n,
            const scalar min,
            const scalar max,
            const scalar mu,
            const scalar sigma
        );

        //- Construct copy
        normal(const normal& d, const label sampleQ);

        //- Construct and return a clone
        virtual autoPtr<distribution> clone(const label sampleQ) const
        {
            return autoPtr<distribution>(new normal(*this, sampleQ));
        }


    //- Destructor
    virtual ~normal();


    // Member Functions

        //- Sample the distribution
        virtual scalar sample() const;

        //- Sample the distribution
        using FieldDistribution<unintegrableForNonZeroQ, normal>::sample;

        //- Return the minimum value
        virtual scalar min() const;

        //- Return the maximum value
        virtual scalar max() const;

        //- Return the mean value
        virtual scalar mean() const;

        //- Return coordinates to plot across the range of the distribution
        virtual tmp<scalarField> x(const label n) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace distributions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
