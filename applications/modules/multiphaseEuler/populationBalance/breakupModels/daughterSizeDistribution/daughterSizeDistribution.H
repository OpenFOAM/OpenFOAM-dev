/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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
    Foam::populationBalance::breakupModels::daughterSizeDistribution

Description
    Base class for breakup models which provide a total breakup rate and a
    separate daughter size distribution function.

See also
    Foam::populationBalanceModel

SourceFiles
    daughterSizeDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef daughterSizeDistribution_H
#define daughterSizeDistribution_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace breakupModels
{

/*---------------------------------------------------------------------------*\
                        Class daughterSizeDistribution Declaration
\*---------------------------------------------------------------------------*/

class daughterSizeDistribution
:
    public breakupModel
{
protected:

    // Protected Data

        //- Daughter size distribution model
        autoPtr<daughterSizeDistributionModel> dsdPtr_;


public:

    // Constructor

        daughterSizeDistribution
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~daughterSizeDistribution()
    {}


    // Member Functions

        //- Return const-reference to daughter size distribution model
        const daughterSizeDistributionModel& dsd() const
        {
            return dsdPtr_();
        }

        //- Return the breakup rate for a group
        virtual tmp<volScalarField::Internal> rate(const label i) const = 0;

        //- Pre-compute diameter independent expressions
        virtual void precompute();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace breakupModels
} // End namespace populationBalance
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
