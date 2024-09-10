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
    Foam::distributions::tabulatedCumulative

Description
    Distribution in which the cumulative distribution function is given as a
    table of values

Usage
    Example usage:
    \verbatim
    {
        type    tabulatedCumulative;
        Q       0;
        distribution
        (
            (10e-6 0)
            (12e-6 0.1972235)
            (15e-6 2.6711166)
            (19e-6 11.643361)
            (25e-6 28.968245)
            (31e-6 52.155796)
            (39e-6 74.381959)
            (50e-6 91.162850)
            (63e-6 100)
        );
    }
    \endverbatim

SourceFiles
    tabulatedCumulative.C

See also
    Foam::distribution

\*---------------------------------------------------------------------------*/

#ifndef tabulatedCumulative_H
#define tabulatedCumulative_H

#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{

/*---------------------------------------------------------------------------*\
                     Class tabulatedCumulative Declaration
\*---------------------------------------------------------------------------*/

class tabulatedCumulative
:
    public FieldDistribution<distribution, tabulatedCumulative>
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
    TypeName("tabulatedCumulative");


    // Constructors

        //- Construct from a dictionary
        tabulatedCumulative
        (
            const unitConversion& units,
            const dictionary& dict,
            const label sampleQ,
            randomGenerator&& rndGen
        );

        //- Construct copy
        tabulatedCumulative(const tabulatedCumulative& d, const label sampleQ);

        //- Construct and return a clone
        virtual autoPtr<distribution> clone(const label sampleQ) const
        {
            return autoPtr<distribution>
            (
                new tabulatedCumulative(*this, sampleQ)
            );
        }


    //- Destructor
    virtual ~tabulatedCumulative();


    // Member Functions

        //- Sample the distribution
        virtual scalar sample() const;

        //- Sample the distribution
        using FieldDistribution<distribution, tabulatedCumulative>::sample;

        //- Return the minimum value
        virtual scalar min() const;

        //- Return the maximum value
        virtual scalar max() const;

        //- Return the mean value
        virtual scalar mean() const;

        //- Return the integral of the PDF multiplied by an integer power of x
        virtual tmp<scalarField> integralPDFxPow
        (
            const scalarField& x,
            const label e,
            const bool consistent = false
        ) const;

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
