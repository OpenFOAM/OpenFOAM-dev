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
    Foam::distributions::unintegrable

Description
    Base class for distributions that do not have a closed integral form for
    the cumulative density function (CDF) for some or all effective size
    exponents.

SourceFiles
    unintegrable.C

See also
    Foam::distribution

\*---------------------------------------------------------------------------*/

#ifndef unintegrable_H
#define unintegrable_H

#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{

/*---------------------------------------------------------------------------*\
                        Class unintegrable Declaration
\*---------------------------------------------------------------------------*/

class unintegrable
:
    public distribution
{
private:

    // Private Data

        //- Number of intervals to use over the range of the distribution
        const label n_;

        //- X-coordinates. These are spaced to contain equal increments of the
        //  CDF. Sampling is a inversion of the CDF; i.e., solve CDF(x) = s for
        //  x, where s is a random sample between 0 and 1. Uniform increments
        //  of CDF make this inversion trivial.
        mutable autoPtr<scalarField> xPtr_;

        //- Values of Phi (i.e., the un-normalised CDF) at the minimum and
        //  maximum x-coordinates
        mutable autoPtr<Pair<scalar>> Phi01Ptr_;

        //- Values of the PDF at the X-coordinates
        mutable autoPtr<scalarField> PDFPtr_;


    // Private Member Functions

        //- X-coordinates
        const scalarField& x() const;

        //- Values of the PDF at the X-coordinates
        const scalarField& PDF() const;


protected:

    // Protected Member Functions

        //- Return values of the un-normalised PDF for the given size exponent
        //  and x-coordinates. Must be provided by derivations.
        virtual tmp<scalarField> phi
        (
            const label q,
            const scalarField& x
        ) const = 0;

        //- Return values of the un-normalised CDF for the given size exponent
        //  and x-coordinates. Can be overloaded by derivations, for example, if
        //  there is a simpler analytic solution for certain effective size
        //  exponents (probably zero).
        virtual tmp<scalarField> Phi
        (
            const label q,
            const scalarField& x
        ) const;

        //- Return values of the un-normalised CDF at the minimum and maximum
        //  x-coordinates for the given size exponent. Can be overloaded by
        //  derivations, for example, if there is a simpler analytic solution
        //  for certain effective size exponents (probably zero).
        virtual Pair<scalar> Phi01(const label q) const;

        //- Access cached values of the un-normalised CDF at the minimum and
        //  maximum x-coordinates.
        const Pair<scalar>& Phi01() const;


public:

    // Static Member Functions

        //- Integrate the values y with respect to the coordinates x
        static tmp<scalarField> integrate
        (
            const scalarField& x,
            const scalarField& y
        );

        //- Integrate the values x*y with respect to the coordinates x
        static tmp<scalarField> integrateX
        (
            const scalarField& x,
            const scalarField& y
        );

        //- Sample an interval, given the interval's bounding x-coordinates,
        //  values Phi (the un-normalised CDF), and a random sample. First
        //  order interpolation.
        static scalar sampleInterval
        (
            const Pair<scalar>& x,
            const Pair<scalar>& Phi,
            const scalar s
        );

        //- Sample an interval, given the interval's bounding x-coordinates,
        //  values phi and Phi (the un-normalised PDF and CDF), and a random
        //  sample. Second order interpolation.
        static scalar sampleInterval
        (
            const Pair<scalar>& x,
            const Pair<scalar>& phi,
            const Pair<scalar>& Phi,
            const scalar s
        );

        //- Sample a discretised distribution, given the x-coordinates,
        //  values Phi (the un-normalised CDF), and a random sample. First
        //  order interpolation.
        static scalar sample
        (
            const scalarField& x,
            const scalarField& Phi,
            const scalar s
        );

        //- Sample a discretised distribution, given the x-coordinates, values
        //  phi and Phi (the un-normalised PDF and CDF), and a random sample.
        //  Second order interpolation.
        static scalar sample
        (
            const scalarField& x,
            const scalarField& phi,
            const scalarField& Phi,
            const scalar s
        );


    // Constructors

        //- Construct from a dictionary
        unintegrable
        (
            const word& name,
            const dictionary& dict,
            randomGenerator& rndGen,
            const label sampleQ
        );

        //- Construct from components
        unintegrable
        (
            randomGenerator& rndGen,
            const label Q,
            const label sampleQ,
            const label n
        );

        //- Construct copy
        unintegrable(const unintegrable& d, const label sampleQ);


    //- Destructor
    virtual ~unintegrable();


    // Member Functions

        //- Sample the distribution
        scalar sample() const;

        //- Sample the distribution
        using distribution::sample;

        //- Return the mean value
        virtual scalar mean() const;

        //- Return coordinates to plot across the range of the distribution
        using distribution::x;

        //- Return the distribution probability density function
        virtual tmp<scalarField> PDF(const scalarField& x) const;
};


/*---------------------------------------------------------------------------*\
                        Class unintegrable Declaration
\*---------------------------------------------------------------------------*/

class unintegrableForNonZeroQ
:
    public unintegrable
{
protected:

    // Protected Member Functions

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace distributions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
