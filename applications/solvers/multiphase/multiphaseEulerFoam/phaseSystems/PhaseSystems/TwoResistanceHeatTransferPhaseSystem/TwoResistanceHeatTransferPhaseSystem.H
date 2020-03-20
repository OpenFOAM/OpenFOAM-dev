/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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
    Foam::TwoResistanceHeatTransferPhaseSystem

Description
    Class which models interfacial heat transfer between a number of phases.
    Two heat transfer models are stored at each interface, one for each phase.
    This permits definition of an interface temperature with which heat transfer
    occurs. It also allows derived systems to define other thermodynamic
    properties at the interface and therefore represent phase changes.

See also
    OneResistanceHeatTransferPhaseSystem

SourceFiles
    TwoResistanceHeatTransferPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef TwoResistanceHeatTransferPhaseSystem_H
#define TwoResistanceHeatTransferPhaseSystem_H

#include "twoResistanceHeatTransferPhaseSystem.H"
#include "HeatTransferPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class modelType>
class BlendedInterfacialModel;

class blendingMethod;
class heatTransferModel;

/*---------------------------------------------------------------------------*\
            Class TwoResistanceHeatTransferPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class TwoResistanceHeatTransferPhaseSystem
:
    public twoResistanceHeatTransferPhaseSystem,
    public HeatTransferPhaseSystem<BasePhaseSystem>
{
protected:

    // Protected typedefs

        typedef HashTable
        <
            Pair<autoPtr<BlendedInterfacialModel<heatTransferModel>>>,
            phasePairKey,
            phasePairKey::hash
        > heatTransferModelTable;

        typedef
            typename HeatTransferPhaseSystem<BasePhaseSystem>::latentHeatScheme
            latentHeatScheme;


    // Protected data

        // Sub Models

            //- Heat transfer models
            heatTransferModelTable heatTransferModels_;


    // Protected Member Functions

        //- Add energy transfer terms which result from bulk mass transfers
        //  and phase changes
        using HeatTransferPhaseSystem<BasePhaseSystem>::addDmdtHefs;

        //- Add energy transfer terms which result from bulk phase changes
        //  that are coupled to the two-resistance heat transfer model
        void addDmdtHefs
        (
            const phaseSystem::dmdtfTable& dmdtfs,
            const phaseSystem::dmdtfTable& Tfs,
            const latentHeatScheme scheme,
            const latentHeatTransfer transfer,
            phaseSystem::heatTransferTable& eqns
        ) const;

        //- Add energy transfer terms which result from specie mass transfers
        //  and phase changes
        using HeatTransferPhaseSystem<BasePhaseSystem>::addDmidtHefs;

        //- Add energy transfer terms which result from specie phase changes
        //  that are coupled to the two-resistance heat transfer model
        void addDmidtHefs
        (
            const phaseSystem::dmidtfTable& dmidtfs,
            const phaseSystem::dmdtfTable& Tfs,
            const latentHeatScheme scheme,
            const latentHeatTransfer transfer,
            phaseSystem::heatTransferTable& eqns
        ) const;


public:

    // Constructors

        //- Construct from fvMesh
        TwoResistanceHeatTransferPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~TwoResistanceHeatTransferPhaseSystem();


    // Member Functions

        //- Return the heat transfer matrices
        virtual autoPtr<phaseSystem::heatTransferTable> heatTransfer() const;

        //- Correct the energy transport e.g. alphat and Tf
        virtual void correctEnergyTransport();

        //- Correct the interface thermodynamics
        virtual void correctInterfaceThermo() = 0;

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TwoResistanceHeatTransferPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
