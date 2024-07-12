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
    Foam::distributions::exponential

Description
    Exponential distribution, scaled so that it spans between a specified
    minimum and maximum value, rather than from zero to infinity

    \f[
        PDF(x) = \lambda \exp(- \lambda x)
    \f]

Usage
    Example usage:
    \verbatim
    {
        type    exponential;
        Q       3;
        min     0.01;
        max     0.5;
        lambda  3;
    }
    \endverbatim

SourceFiles
    exponential.C

See also
    Foam::distribution

\*---------------------------------------------------------------------------*/

#ifndef exponential_H
#define exponential_H

#include "unintegrableForNonZeroQ.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{

/*---------------------------------------------------------------------------*\
                        Class exponential Declaration
\*---------------------------------------------------------------------------*/

class exponential
:
    public FieldDistribution<unintegrableForNonZeroQ, exponential>
{
    // Private Data

        //- Minimum value
        const scalar min_;

        //- Maximum value
        const scalar max_;

        //- Parameter
        const scalar lambda_;


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
    TypeName("exponential");


    // Constructors

        //- Construct from a dictionary
        exponential
        (
            const unitConversion& units,
            const dictionary& dict,
            const label sampleQ,
            randomGenerator&& rndGen
        );

        //- Construct copy
        exponential(const exponential& d, const label sampleQ);

        //- Construct and return a clone
        virtual autoPtr<distribution> clone(const label sampleQ) const
        {
            return autoPtr<distribution>(new exponential(*this, sampleQ));
        }


    //- Destructor
    virtual ~exponential();


    // Member Functions

        //- Sample the distribution for zero effective size exponent
        virtual scalar sampleForZeroQ() const;

        //- Sample the distribution
        using FieldDistribution<unintegrableForNonZeroQ, exponential>::sample;

        //- Return the minimum value
        virtual scalar min() const;

        //- Return the maximum value
        virtual scalar max() const;

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
