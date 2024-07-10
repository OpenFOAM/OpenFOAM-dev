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
    Foam::distributions::uniform

Description
    Distribution in which all values between a specified minimum and maximum
    have the same probability

Usage
    Example usage:
    \verbatim
    {
        type    uniform;
        Q       0;
        min     0.5;
        max     1;
    }
    \endverbatim

SourceFiles
    uniform.C

See also
    Foam::distribution

\*---------------------------------------------------------------------------*/

#ifndef uniform_H
#define uniform_H

#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{

/*---------------------------------------------------------------------------*\
                           Class uniform Declaration
\*---------------------------------------------------------------------------*/

class uniform
:
    public FieldDistribution<distribution, uniform>
{
    // Private Data

        //- Distribution minimum
        const scalar min_;

        //- Distribution maximum
        const scalar max_;

        //- Cached indefinite integral of the PDF at the minimum value
        const scalar Phi0_;

        //- Cached indefinite integral of the PDF at the maximum value
        const scalar Phi1_;


    // Private Member Functions

        //- Calculate the indefinite integral of the PDF at the given values
        template<class Type>
        auto Phi(const Type& x) const;


public:

    //- Runtime type information
    TypeName("uniform");


    // Constructors

        //- Construct from a dictionary
        uniform
        (
            const unitConversion& units,
            const dictionary& dict,
            const label sampleQ,
            randomGenerator&& rndGen
        );

        //- Construct copy
        uniform(const uniform& d, const label sampleQ);

        //- Construct and return a clone
        virtual autoPtr<distribution> clone(const label sampleQ) const
        {
            return autoPtr<distribution>(new uniform(*this, sampleQ));
        }


    //- Destructor
    virtual ~uniform();


    // Member Functions

        //- Sample the distribution
        virtual scalar sample() const;

        //- Sample the distribution
        using FieldDistribution<distribution, uniform>::sample;

        //- Return the minimum value
        virtual scalar min() const;

        //- Return the maximum value
        virtual scalar max() const;

        //- Return the mean value
        virtual scalar mean() const;

        //- Return the cumulative density function at the given coordinates
        virtual tmp<scalarField> CDF(const scalarField& x) const;

        //- Write to a stream
        virtual void write(Ostream& os, const unitConversion& units) const;

        //- Return values to plot the probability density function
        virtual tmp<scalarField> plotPDF(const scalarField& x) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace distributions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
