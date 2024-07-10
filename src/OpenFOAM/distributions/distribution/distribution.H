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
#include "randomGenerator.H"
#include "fieldTypes.H"
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

        //- Distribution size exponent
        const label Q_;

        //- Sample size exponent
        const label sampleQ_;

        //- Random number generator
        mutable randomGenerator rndGen_;


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

        //- Sample the distribution into components of a primitive type
        #define VIRTUAL_SAMPLE_TYPE(Type, nullArg)                             \
            virtual Type CAT(sample, CAPITALIZE(Type))() const = 0;
        FOR_ALL_FIELD_TYPES(VIRTUAL_SAMPLE_TYPE);
        #undef VIRTUAL_SAMPLE_TYPE


public:

    //- Runtime type information
    TypeName("distribution");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        distribution,
        dictionary,
        (
            const unitConversion& units,
            const dictionary& dict,
            const label sampleQ,
            randomGenerator&& rndGen
        ),
        (units, dict, sampleQ, std::move(rndGen))
    );


    // Constructors

        //- Construct from dictionary
        distribution
        (
            const word& name,
            const unitConversion& units,
            const dictionary& dict,
            const label sampleQ,
            randomGenerator&& rndGen
        );

        //- Construct from components
        distribution
        (
            const label Q,
            const label sampleQ,
            randomGenerator&& rndGen
        );

        //- Construct copy
        distribution(const distribution& d, const label sampleQ);

        //- Construct and return a clone
        virtual autoPtr<distribution> clone(const label sampleQ) const = 0;

        //- Construct and return a clone
        inline autoPtr<distribution> clone() const
        {
            return clone(sampleQ_);
        }


    // Selectors

        //- Select from dictionary and a random generator
        static autoPtr<distribution> New
        (
            const unitConversion& units,
            const dictionary& dict,
            const label sampleQ,
            randomGenerator&& rndGen,
            const bool report = true
        );

        //- Select from a dictionary and a random generator seed and global flag
        static autoPtr<distribution> New
        (
            const unitConversion& units,
            const dictionary& dict,
            const label sampleQ,
            const randomGenerator::seed& s,
            const bool global = false,
            const bool report = true
        );

        //- Re-select with a different sample size exponent
        static autoPtr<distribution> New
        (
            autoPtr<distribution>& dPtr,
            const label sampleQ
        );


    //- Destructor
    virtual ~distribution();


    // Member Functions

        //- Sample the distribution
        virtual scalar sample() const = 0;

        //- Sample the distribution into components of a primitive type
        template<class Type>
        Type sample() const;

        //- Sample the distribution into a field
        virtual tmp<scalarField> sample(const label n) const = 0;

        //- Return the minimum value
        virtual scalar min() const = 0;

        //- Return the maximum value
        virtual scalar max() const = 0;

        //- Return the mean value
        virtual scalar mean() const = 0;

        //- Return the cumulative density function at the given coordinates
        virtual tmp<scalarField> CDF(const scalarField& x) const = 0;

        //- Write to a stream
        virtual void write(Ostream& os, const unitConversion& units) const;

        //- Write the state to a stream
        virtual void writeState(Ostream& os) const;

        //- Return coordinates to plot across the range of the distribution
        virtual tmp<scalarField> plotX(const label n) const;

        //- Return values to plot the probability density function
        virtual tmp<scalarField> plotPDF(const scalarField& x) const = 0;
};


#define DISTRIBUTION_TEMPLATED_SAMPLE_TYPE(Type, nullArg)                      \
    template<>                                                                 \
    inline Type Foam::distribution::sample<Type>() const                       \
    {                                                                          \
        return CAT(sample, CAPITALIZE(Type))();                                \
    }
FOR_ALL_FIELD_TYPES(DISTRIBUTION_TEMPLATED_SAMPLE_TYPE);
#undef DISTRIBUTION_TEMPLATED_SAMPLE_TYPE


void writeEntry
(
    Ostream& os,
    const unitConversion&,
    const distribution& d,
    const bool write = true,
    const bool writeState = true
);


/*---------------------------------------------------------------------------*\
                      Class FieldDistribution Declaration
\*---------------------------------------------------------------------------*/

template<class Base, class Derived>
class FieldDistribution
:
    public Base
{
protected:

    // Protected Member Functions

        //- Sample the distribution into components of a primitive type
        #define VIRTUAL_SAMPLE_TYPE(Type, nullArg)                             \
            virtual Type CAT(sample, CAPITALIZE(Type))() const                 \
            {                                                                  \
                return sample<Type>();                                         \
            }
        FOR_ALL_FIELD_TYPES(VIRTUAL_SAMPLE_TYPE);
        #undef VIRTUAL_SAMPLE_TYPE


public:

    // Constructors

        //- Inherit constructors
        using Base::Base;


    // Member Functions

        //- Sample the distribution
        using Base::sample;

        //- Sample the distribution into components of a primitive type
        template<class Type>
        Type sample() const;

        //- Sample the distribution into a field
        virtual tmp<scalarField> sample(const label n) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "distributionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
