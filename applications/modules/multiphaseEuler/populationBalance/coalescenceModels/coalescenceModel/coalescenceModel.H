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
    Foam::populationBalance::coalescenceModel

Description
    Base class for coalescence models.

See also
    Foam::populationBalanceModel

SourceFiles
    coalescenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef coalescenceModel_H
#define coalescenceModel_H

#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{

/*---------------------------------------------------------------------------*\
                      Class coalescenceModel Declaration
\*---------------------------------------------------------------------------*/

class coalescenceModel
{
protected:

    // Protected Data

        //- Reference to the populationBalanceModel
        const populationBalanceModel& popBal_;


public:

    //- Runtime type information
    TypeName("coalescenceModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            coalescenceModel,
            dictionary,
            (
                const populationBalanceModel& popBal,
                const dictionary& dict
            ),
            (popBal, dict)
        );


    // Constructor

        //- Construct from a population balance and its dictionary
        coalescenceModel
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Selector
    static autoPtr<coalescenceModel> New
    (
        const populationBalanceModel& popBal,
        const dictionary& dict
    );


    //- Destructor
    virtual ~coalescenceModel()
    {}


    // Member Functions

        //- Does this model coalesce?
        virtual bool coalesces() const;

        //- Pre-compute diameter independent expressions
        virtual void precompute();

        //- Return the coalescence rate between two groups
        virtual tmp<volScalarField::Internal> rate
        (
            const label i,
            const label j
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalance
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
