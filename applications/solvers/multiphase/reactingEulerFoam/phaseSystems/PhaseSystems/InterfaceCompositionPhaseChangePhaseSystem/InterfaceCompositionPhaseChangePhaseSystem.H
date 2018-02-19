/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::InterfaceCompositionPhaseChangePhaseSystem

Description
    Class to provide interfacial heat and mass transfer between a number of
    phases according to a interface composition model.

    The interface temperature is calculated such that the net rate at which the
    heat is transferred to the interface is equal to the latent heat consumed by
    the mass transfer.

SourceFiles
    InterfaceCompositionPhaseChangePhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef InterfaceCompositionPhaseChangePhaseSystem_H
#define InterfaceCompositionPhaseChangePhaseSystem_H

#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class interfaceCompositionModel;

/*---------------------------------------------------------------------------*\
                 Class InterfaceCompositionPhaseChangePhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class InterfaceCompositionPhaseChangePhaseSystem
:
    public BasePhaseSystem
{
protected:

    // Protected typedefs

        typedef HashTable
        <
            autoPtr<interfaceCompositionModel>,
            phasePairKey,
            phasePairKey::hash
        > interfaceCompositionModelTable;


    // Protected data

        // Sub Models

            //- Interface composition models
            interfaceCompositionModelTable interfaceCompositionModels_;

            //- Interfacial Mass transfer rate
            HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>
                iDmdt_;

            //- Explicit part of the mass transfer rate
            HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>
                iDmdtExplicit_;

    // Protected member functions

        //- Add the mass-based momentum transfer to the equations
        void addMomentumTransfer
        (
            phaseSystem::momentumTransferTable& eqns
        ) const;


public:

    // Constructors

        //- Construct from fvMesh
        InterfaceCompositionPhaseChangePhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~InterfaceCompositionPhaseChangePhaseSystem();


    // Member Functions

        //- Return the interfacial mass flow rate
        virtual tmp<volScalarField> iDmdt(const phasePairKey& key) const;

        //- Return the total interfacial mass transfer rate for phase
        virtual tmp<volScalarField> iDmdt(const phaseModel& phase) const;

        //- Return the interfacial mass flow rate
        virtual tmp<volScalarField> dmdt(const phasePairKey& key) const
        {
            return BasePhaseSystem::dmdt(key) + this->iDmdt(key);
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

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "InterfaceCompositionPhaseChangePhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
