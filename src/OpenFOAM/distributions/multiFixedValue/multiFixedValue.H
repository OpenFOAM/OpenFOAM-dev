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
    Foam::distributions::multiFixedValue

Description
    Distribution which comprises a list of fixed values with given
    probabilities. The probabilities are normalised so that they sum to one, so
    they can be provided as fractions or percentages or with any other
    arbitrary scale factor.

Usage
    Example usage:
    \verbatim
    {
        type    multiFixedValue;
        values
        (
            (1e-6 30)
            (2e-6 60)
            (5e-6 10)
        );
    }
    \endverbatim

SourceFiles
    multiFixedValue.C

See also
    Foam::distribution

\*---------------------------------------------------------------------------*/

#ifndef multiFixedValue_H
#define multiFixedValue_H

#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{

/*---------------------------------------------------------------------------*\
                         Class multiFixedValue Declaration
\*---------------------------------------------------------------------------*/

class multiFixedValue
:
    public FieldDistribution<distribution, multiFixedValue>
{
    // Private Data

        //- X-coordinates
        scalarField x_;

        //- Probabilities
        scalarField P_;

        //- Cumulative sum of the probabilities
        scalarField sumP_;


public:

    //- Runtime type information
    TypeName("multiFixedValue");


    // Constructors

        //- Construct from a dictionary
        multiFixedValue
        (
            const unitConversion& units,
            const dictionary& dict,
            const label sampleQ,
            randomGenerator&& rndGen
        );

        //- Construct copy
        multiFixedValue(const multiFixedValue& d, const label sampleQ);

        //- Construct and return a clone
        virtual autoPtr<distribution> clone(const label sampleQ) const
        {
            return autoPtr<distribution>(new multiFixedValue(*this, sampleQ));
        }


    //- Destructor
    virtual ~multiFixedValue();


    // Member Functions

        //- Sample the distribution
        virtual scalar sample() const;

        //- Sample the distribution
        using FieldDistribution<distribution, multiFixedValue>::sample;

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

        //- Return coordinates to plot across the range of the distribution
        virtual tmp<scalarField> plotX(const label n) const;

        //- Return values to plot the probability density function
        virtual tmp<scalarField> plotPDF(const scalarField& x) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace distributions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
