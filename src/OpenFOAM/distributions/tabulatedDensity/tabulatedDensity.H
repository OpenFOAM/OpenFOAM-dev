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
    Foam::distributions::tabulatedDensity

Description
    Distribution in which the probability density function is given as a table
    of values

Usage
    Example usage:
    \verbatim
    {
        type    tabulatedDensity;
        Q       0;
        distribution
        (
            (10e-6 0.0025)
            (12e-6 0.2795)
            (15e-6 2.3988)
            (19e-6 6.3888)
            (25e-6 10.3153)
            (31e-6 12.0030)
            (39e-6 10.8427)
            (50e-6 6.1333)
            (63e-6 3.4688)
        );

    }
    \endverbatim

SourceFiles
    tabulatedDensity.C

See also
    Foam::distribution

\*---------------------------------------------------------------------------*/

#ifndef tabulatedDensity_H
#define tabulatedDensity_H

#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{

/*---------------------------------------------------------------------------*\
                        Class tabulatedDensity Declaration
\*---------------------------------------------------------------------------*/

class tabulatedDensity
:
    public FieldDistribution<distribution, tabulatedDensity>
{
private:

    // Private Data

        //- X-coordinates
        scalarField x_;

        //- Probability Density Function
        scalarField PDF_;

        //- Cumulative Distribution Function
        scalarField CDF_;


public:

    //- Runtime type information
    TypeName("tabulatedDensity");


    // Constructors

        //- Construct from a dictionary
        tabulatedDensity
        (
            const unitConversion& units,
            const dictionary& dict,
            const label sampleQ,
            randomGenerator&& rndGen
        );

        //- Construct copy
        tabulatedDensity(const tabulatedDensity& d, const label sampleQ);

        //- Construct and return a clone
        virtual autoPtr<distribution> clone(const label sampleQ) const
        {
            return autoPtr<distribution>(new tabulatedDensity(*this, sampleQ));
        }


    //- Destructor
    virtual ~tabulatedDensity();


    // Member Functions

        //- Sample the distribution
        virtual scalar sample() const;

        //- Sample the distribution
        using FieldDistribution<distribution, tabulatedDensity>::sample;

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
