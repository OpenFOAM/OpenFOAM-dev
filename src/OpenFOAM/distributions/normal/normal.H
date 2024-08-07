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
    Foam::distributions::normal

Description
    Normal distribution, scaled so that it spans between a specified
    minimum and maximum value, rather than from zero to infinity

    \f[
        PDF(x) = \frac{1}{\sigma \sqrt{2 \pi}} \exp \\
        \left( - \frac{1}{2} \left( \frac{x - \mu}{\sigma} \right)^2 \right)
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

#include "unintegrableForNonZeroQ.H"

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


    // Private Member Functions

        //- Return values of the un-normalised PDF for the given size exponent
        //  and x-coordinates.
        virtual tmp<scalarField> phi
        (
            const label q,
            const scalarField& x
        ) const;

        //- Return values of the un-normalised CDF for zero effective size
        //  exponent and given x-coordinates
        virtual tmp<scalarField> PhiForZeroQ(const scalarField& x) const;


public:

    //- Runtime type information
    TypeName("normal");


    //- Permit the multiNormal distribution to use private parts of this class
    friend class multiNormal;


    // Constructors

        //- Construct from a dictionary
        normal
        (
            const unitConversion& units,
            const dictionary& dict,
            const label sampleQ,
            randomGenerator&& rndGen
        );

        //- Construct from components
        normal
        (
            const label Q,
            const label sampleQ,
            randomGenerator&& rndGen,
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

        //- Sample the distribution for zero effective size exponent
        virtual scalar sampleForZeroQ() const;

        //- Sample the distribution for zero effective size exponent
        scalar sampleForZeroQ(const scalar s) const;

        //- Sample the distribution
        using FieldDistribution<unintegrableForNonZeroQ, normal>::sample;

        //- Return the minimum value
        virtual scalar min() const;

        //- Return the maximum value
        virtual scalar max() const;

        //- Return the mean value
        virtual scalar mean() const;

        //- Return the mean value
        inline scalar mu() const
        {
            return mu_;
        }

        //- Return the standard deviation
        inline scalar sigma() const
        {
            return sigma_;
        }

        //- Write to a stream
        virtual void write(Ostream& os, const unitConversion& units) const;

        //- Return coordinates to plot across the range of the distribution
        virtual tmp<scalarField> plotX(const label n) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace distributions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
