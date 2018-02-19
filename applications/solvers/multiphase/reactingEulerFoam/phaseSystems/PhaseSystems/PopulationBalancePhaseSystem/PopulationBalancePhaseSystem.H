/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::PopulationBalancePhaseSystem

Description
    Class which provides population balance functionality.

See also
    Foam::diameterModels::populationBalanceModel

SourceFiles
    PopulationBalancePhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef PopulationBalancePhaseSystem_H
#define PopulationBalancePhaseSystem_H

#include "HeatAndMassTransferPhaseSystem.H"
#include "saturationModel.H"
#include "Switch.H"

#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class PopulationBalancePhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class PopulationBalancePhaseSystem
:
    public BasePhaseSystem
{

protected:

    // Protected data

        //- populationBalanceModels
        PtrList<diameterModels::populationBalanceModel> populationBalances_;

        //- Interfacial Mass transfer rate
        HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash> pDmdt_;


    // Protected member functions

        //- Add the mass-based momentum transfer to the equations
        void addMomentumTransfer
        (
            phaseSystem::momentumTransferTable& eqns
        ) const;


public:

    // Constructors

        //- Construct from fvMesh
        PopulationBalancePhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~PopulationBalancePhaseSystem();


    // Member Functions

        //- Return the interfacial mass flow rate
        virtual tmp<volScalarField> pDmdt(const phasePairKey& key) const;

        //- Return the interfacial mass flow rate
        virtual tmp<volScalarField> dmdt(const phasePairKey& key) const
        {
            return BasePhaseSystem::dmdt(key) + this->pDmdt(key);
        };

        //- Return the total interfacial mass transfer rate for phase
        virtual tmp<volScalarField> dmdt(const phaseModel& phase) const;

        //- Return the momentum transfer matrices
        virtual autoPtr<phaseSystem::momentumTransferTable>
            momentumTransfer() const;

        //- Return the momentum transfer matrices for the face based algorithm
        virtual autoPtr<phaseSystem::momentumTransferTable>
            momentumTransferf() const;

        //- Return the heat transfer matrices
        virtual autoPtr<phaseSystem::heatTransferTable> heatTransfer() const;

        //- Return the mass transfer matrices
        virtual autoPtr<phaseSystem::massTransferTable> massTransfer() const;

        //- Read base phaseProperties dictionary
        virtual bool read();

        //- Solve all population balance equations
        virtual void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PopulationBalancePhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
