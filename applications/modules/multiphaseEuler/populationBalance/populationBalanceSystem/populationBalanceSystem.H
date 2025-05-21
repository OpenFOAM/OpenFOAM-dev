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
    Foam::populationBalanceSystem

Description
    Class which provides population balance functionality. Holds a number of
    population balances and provides access to the transfers that result from
    the population balance interactions.

See also
    Foam::diameterModels::populationBalanceModel

SourceFiles
    populationBalanceSystem.C

\*---------------------------------------------------------------------------*/

#ifndef populationBalanceSystem_H
#define populationBalanceSystem_H

#include "phaseSystem.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class populationBalanceSystem Declaration
\*---------------------------------------------------------------------------*/

class populationBalanceSystem
{
    // Private Data

        //- Reference to the phase system
        const phaseSystem& fluid_;

        //- Population balances
        PtrList<diameterModels::populationBalanceModel> populationBalances_;


    // Private Member Functions

        //- Add interface mass transfer rates to phase-mass transfer rates
        void addDmdts
        (
            const diameterModels::populationBalanceModel::dmdtfTable& dmdtfs,
            PtrList<volScalarField::Internal>& dmdts
        ) const;

        //- Add momentum transfer terms which result from mass transfers
        void addDmdtUfs
        (
            const diameterModels::populationBalanceModel::dmdtfTable& dmdtfs,
            HashPtrTable<fvVectorMatrix>& eqns
        ) const;

        //- Add energy transfer terms which result from mass transfers
        void addDmdtHefs
        (
            const diameterModels::populationBalanceModel::dmdtfTable& dmdtfs,
            HashPtrTable<fvScalarMatrix>& eqns
        ) const;

        //- Add specie transfer terms which result from mass transfers
        void addDmdtYfs
        (
            const diameterModels::populationBalanceModel::dmdtfTable& dmdtfs,
            HashPtrTable<fvScalarMatrix>& eqns
        ) const;


public:

    //- Runtime type information
    TypeName("populationBalanceSystem");


    // Constructors

        //- Construct from a phase system
        populationBalanceSystem(const phaseSystem&);


    //- Destructor
    virtual ~populationBalanceSystem();


    // Member Functions

        //- Return the mass transfer rates for each phase
        PtrList<volScalarField::Internal> dmdts() const;

        //- Return the momentum transfer matrices for the cell-based algorithm
        autoPtr<HashPtrTable<fvVectorMatrix>> momentumTransfer();

        //- Return the momentum transfer matrices for the face-based algorithm
        autoPtr<HashPtrTable<fvVectorMatrix>> momentumTransferf();

        //- Return the heat transfer matrices
        autoPtr<HashPtrTable<fvScalarMatrix>> heatTransfer() const;

        //- Return the specie transfer matrices
        autoPtr<HashPtrTable<fvScalarMatrix>> specieTransfer() const;

        //- Solve all population balance equations
        void solve();

        //- Correct the population balances
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
