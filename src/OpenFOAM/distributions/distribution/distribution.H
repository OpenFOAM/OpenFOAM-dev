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
    Foam::distribution

Description
    Base class for statistical distributions

    All distributions (except fixedValue) require a "size exponent", Q, to be
    specified along with their other coefficients. If a distribution's CDF(x)
    (cumulative distribution function) represents what proportion of the
    distribution takes a value below x, then Q determines what is meant by
    "proportion":

    - If Q=0, then "proportion" means the number of sampled values expected to
      be below x divided by the total number of sampled values.

    - If Q=3, then "proportion" means the expected sum of sampled values cubed
      for values below x divided by the total sum of values cubed. If x is a
      length, then this can be interpreted as a proportion of the total volume
      of sampled objects.

    - If Q=2, and x is a length, then the distribution might represent the
      proportion of surface area, and so on...

    In addition to the user-specification of Q defining what size the given
    distribution relates to, an implementation that uses a distribution can
    also programmatically define a samplingQ to determine what sort of sample
    is being constructed; whether the samples should have an equal number
    (sampleQ=0), volume (sampleQ=3), area (sampleQ=2), etc...

SourceFiles
    distribution.C
    distributionNew.C

\*---------------------------------------------------------------------------*/

#ifndef distribution_H
#define distribution_H

#include "dictionary.H"
#include "Random.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class distribution Declaration
\*---------------------------------------------------------------------------*/

class distribution
{
protected:

    // Protected data

        //- Reference to a random number generator
        Random& rndGen_;

        //- Distribution size exponent
        const label Q_;

        //- Sample size exponent
        const label sampleQ_;


    // Protected Member Functions

        //- Validate that the bounds are monotonic
        virtual void validateBounds(const dictionary& dict) const;

        //- Validate that the lower bound is positive
        virtual void validatePositive(const dictionary& dict) const;

        //- Clip the PDF values to zero outside the bounds
        tmp<scalarField> clipPDF
        (
            const scalarField& x,
            const tmp<scalarField>& pdf
        ) const;

        //- Return the effective distribution size exponent
        inline label q() const
        {
            return sampleQ_ - Q_;
        }


public:

    //-Runtime type information
    TypeName("distribution");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        distribution,
        dictionary,
        (
            const dictionary& dict,
            Random& rndGen,
            const label sampleQ
        ),
        (dict, rndGen, sampleQ)
    );


    // Constructors

        //- Construct from dictionary
        distribution
        (
            const word& name,
            const dictionary& dict,
            Random& rndGen,
            const label sampleQ
        );

        //- Construct from components
        distribution(Random& rndGen, const label Q, const label sampleQ);

        //- Construct copy
        distribution(const distribution& d, const label sampleQ);

        //- Construct and return a clone
        virtual autoPtr<distribution> clone(const label sampleQ) const = 0;

        //- Construct and return a clone
        inline autoPtr<distribution> clone() const
        {
            return clone(sampleQ_);
        }


    //- Selector
    static autoPtr<distribution> New
    (
        const dictionary& dict,
        Random& rndGen,
        const label sampleQ
    );


    //- Destructor
    virtual ~distribution();


    // Member Functions

        //- Sample the distribution
        virtual scalar sample() const = 0;

        //- Sample the distribution
        virtual tmp<scalarField> sample(const label n) const = 0;

        //- Return the minimum value
        virtual scalar min() const = 0;

        //- Return the maximum value
        virtual scalar max() const = 0;

        //- Return the mean value
        virtual scalar mean() const = 0;

        //- Report
        void report() const;

        //- Return coordinates to plot across the range of the distribution
        virtual tmp<scalarField> x(const label n) const;

        //- Return the distribution probability density function
        virtual tmp<scalarField> PDF(const scalarField& x) const = 0;
};


/*---------------------------------------------------------------------------*\
                      Class FieldDistribution Declaration
\*---------------------------------------------------------------------------*/

template<class Base, class Derived>
class FieldDistribution
:
    public Base
{
public:

    // Constructors

        //- Inherit constructors
        using Base::Base;


    // Member Functions

        //- Sample the distribution
        virtual tmp<scalarField> sample(const label n) const
        {
            tmp<scalarField> tResult(new scalarField(n));
            scalarField& result = tResult.ref();

            forAll(result, i)
            {
                result[i] =
                    static_cast<const Derived&>(*this).sample();
            }

            return tResult;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
