/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::distributionModel

Description
    A library of runtime-selectable distribution models.

    Returns a sampled value given the expectation (nu) and variance (sigma^2)

    Current distribution models include:
    - exponential
    - fixedValue
    - general
    - multi-normal
    - normal
    - Rosin-Rammler
    - uniform

    The distributionModel is tabulated in equidistant nPoints, in an interval.
    These values are integrated to obtain the cumulated distribution model,
    which is then used to change the distribution from unifrom to
    the actual distributionModel.

SourceFiles
    distributionModel.C
    distributionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef distributionModel_H
#define distributionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "cachedRandom.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{

/*---------------------------------------------------------------------------*\
                     Class distributionModel Declaration
\*---------------------------------------------------------------------------*/

class distributionModel
{

protected:

    // Protected data

        //- Coefficients dictionary
        const dictionary distributionModelDict_;

        //- Reference to the random number generator
        cachedRandom& rndGen_;


    // Protected Member Functions

        //- Check that the distribution model is valid
        virtual void check() const;


public:

    //-Runtime type information
    TypeName("distributionModel");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        distributionModel,
        dictionary,
        (
            const dictionary& dict,
            cachedRandom& rndGen
        ),
        (dict, rndGen)
    );


    // Constructors

        //- Construct from dictionary
        distributionModel
        (
            const word& name,
            const dictionary& dict,
            cachedRandom& rndGen
        );

        //- Construct copy
        distributionModel(const distributionModel& p);

        //- Construct and return a clone
        virtual autoPtr<distributionModel> clone() const
        {
            return autoPtr<distributionModel>(new distributionModel(*this));
        }


    //- Selector
    static autoPtr<distributionModel> New
    (
        const dictionary& dict,
        cachedRandom& rndGen
    );


    //- Destructor
    virtual ~distributionModel();


    // Member Functions

        //- Sample the distributionModel
        virtual scalar sample() const;

        //- Return the minimum value
        virtual scalar minValue() const;

        //- Return the maximum value
        virtual scalar maxValue() const;

        //- Return the maximum value
        virtual scalar meanValue() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace distributionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
