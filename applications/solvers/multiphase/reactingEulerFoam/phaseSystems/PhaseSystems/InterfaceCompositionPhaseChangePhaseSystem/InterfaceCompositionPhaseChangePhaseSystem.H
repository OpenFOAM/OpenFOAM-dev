/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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
class diffusiveMassTransferModel;

/*---------------------------------------------------------------------------*\
                 Class InterfaceCompositionPhaseChangePhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class InterfaceCompositionPhaseChangePhaseSystem
:
    public BasePhaseSystem
{
    // Private typedefs

        typedef HashTable
        <
            Pair<autoPtr<interfaceCompositionModel>>,
            phasePairKey,
            phasePairKey::hash
        > interfaceCompositionModelTable;

        typedef HashTable
        <
            Pair<autoPtr<BlendedInterfacialModel<diffusiveMassTransferModel>>>,
            phasePairKey,
            phasePairKey::hash
        > diffusiveMassTransferModelTable;


    // Private data

        // Sub Models

            //- The number of interface correctors
            const label nInterfaceCorrectors_;

            //- Mass transfer models
            diffusiveMassTransferModelTable diffusiveMassTransferModels_;

            //- Interface composition models
            interfaceCompositionModelTable interfaceCompositionModels_;

            //- The explicit part of the interfacial mass transfer rates
            phaseSystem::dmidtfTable dmidtfSus_;

            //- The implicit part of the interfacial mass transfer rates
            phaseSystem::dmidtfTable dmidtfSps_;


    // Private member functions

        //- Return mass transfers across each interface
        autoPtr<phaseSystem::dmdtfTable> dmdtfs() const;

        //- Return species mass transfers across each interface
        autoPtr<phaseSystem::dmidtfTable> dmidtfs() const;


public:

    // Constructors

        //- Construct from fvMesh
        InterfaceCompositionPhaseChangePhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~InterfaceCompositionPhaseChangePhaseSystem();


    // Member Functions

        //- Return the mass transfer rate for an interface
        virtual tmp<volScalarField> dmdtf(const phasePairKey& key) const;

        //- Return the mass transfer rates for each phase
        virtual PtrList<volScalarField> dmdts() const;

        //- Return the momentum transfer matrices for the cell-based algorithm
        virtual autoPtr<phaseSystem::momentumTransferTable> momentumTransfer();

        //- Return the momentum transfer matrices for the face-based algorithm
        virtual autoPtr<phaseSystem::momentumTransferTable> momentumTransferf();

        //- Return the heat transfer matrices
        virtual autoPtr<phaseSystem::heatTransferTable> heatTransfer() const;

        //- Return the specie transfer matrices
        virtual autoPtr<phaseSystem::specieTransferTable>
            specieTransfer() const;

        //- Correct the fluid properties other than those listed below
        virtual void correct();

        //- Correct the interface temperatures
        virtual void correctInterfaceThermo();

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
