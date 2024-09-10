/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::distributions::standardNormal

Description
    Standard normal distribution. Not selectable.

    \f[
        PDF(x) = \frac{1}{\sqrt{2 \pi}} \exp \left( - \frac{1}{2} x^2 \right)
    \f]

SourceFiles
    standardNormal.C

See also
    Foam::distribution

\*---------------------------------------------------------------------------*/

#ifndef standardNormal_H
#define standardNormal_H

#include "unintegrable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{

/*---------------------------------------------------------------------------*\
                        Class standardNormal Declaration
\*---------------------------------------------------------------------------*/

class standardNormal
:
    public FieldDistribution<distribution, standardNormal>
{
    // Private Data

        //- Constant for approximate error function and inverse error function
        static const scalar a_;


public:

    //- Runtime type information
    TypeName("standardNormal");


    //- Permit the multiNormal distribution to use private parts of this class
    friend class multiNormal;


    // Static Member Functions

        //- Approximate error function
        static tmp<scalarField> approxErf(const scalarField& x);

        //- Approximate error function inverse
        static scalar approxErfInv(const scalar y);


    // Constructors

        //- Construct from a random generator
        standardNormal(randomGenerator&& rndGen);

        //- Construct from a seed
        standardNormal
        (
            const randomGenerator::seed& s,
            const bool global = false
        );

        //- Construct copy
        standardNormal(const standardNormal& d);

        //- Construct and return a clone
        virtual autoPtr<distribution> clone(const label sampleQ) const
        {
            return autoPtr<distribution>(new standardNormal(*this));
        }


    //- Destructor
    virtual ~standardNormal();


    // Member Functions

        //- Sample the distribution
        virtual scalar sample() const;

        //- Sample the distribution
        using FieldDistribution<distribution, standardNormal>::sample;

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
