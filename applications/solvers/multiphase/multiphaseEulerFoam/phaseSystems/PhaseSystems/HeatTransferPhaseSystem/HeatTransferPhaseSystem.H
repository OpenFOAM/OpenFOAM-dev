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
    Foam::HeatTransferPhaseSystem

Description
    ...

SourceFiles
    HeatTransferPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef HeatTransferPhaseSystem_H
#define HeatTransferPhaseSystem_H

#include "heatTransferPhaseSystem.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class HeatTransferPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class HeatTransferPhaseSystem
:
    public heatTransferPhaseSystem,
    public BasePhaseSystem
{
protected:

    // Protected Member Functions

        //- Add energy transfer terms which result from bulk mass transfers
        void addDmdtHefs
        (
            const phaseSystem::dmdtfTable& dmdtfs,
            phaseSystem::heatTransferTable& eqns
        ) const;

        //- Add energy transfer terms which result from specie mass transfers
        void addDmidtHefs
        (
            const phaseSystem::dmidtfTable& dmidtfs,
            phaseSystem::heatTransferTable& eqns
        ) const;

        //- Add energy transfer terms which result from bulk phase changes,
        //  without the latent heat contribution
        void addDmdtHefsWithoutL
        (
            const phaseSystem::dmdtfTable& dmdtfs,
            const phaseSystem::dmdtfTable& Tfs,
            const latentHeatScheme scheme,
            phaseSystem::heatTransferTable& eqns
        ) const;

        //- Add latent heat terms which result from bulk phase changes.
        //  Weight is the proportion of the latent heat contribution applied to
        //  the downwind side of the mass transfer.
        void addDmdtL
        (
            const phaseSystem::dmdtfTable& dmdtfs,
            const phaseSystem::dmdtfTable& Tfs,
            const scalar weight,
            const latentHeatScheme scheme,
            phaseSystem::heatTransferTable& eqns
        ) const;

        //- Add energy transfer terms which result from bulk phase changes.
        //  Weight is the proportion of the latent heat contribution applied to
        //  the downwind side of the mass transfer.
        void addDmdtHefs
        (
            const phaseSystem::dmdtfTable& dmdtfs,
            const phaseSystem::dmdtfTable& Tfs,
            const scalar weight,
            const latentHeatScheme scheme,
            phaseSystem::heatTransferTable& eqns
        ) const;

        //- Add energy transfer terms which result from specie phase changes,
        //  without the latent heat contribution
        void addDmidtHefsWithoutL
        (
            const phaseSystem::dmidtfTable& dmidtfs,
            const phaseSystem::dmdtfTable& Tfs,
            const latentHeatScheme scheme,
            phaseSystem::heatTransferTable& eqns
        ) const;

        //- Add latent heat terms which result from specie phase changes.
        //  Weight is the proportion of the latent heat contribution applied to
        //  the downwind side of the mass transfer.
        void addDmidtL
        (
            const phaseSystem::dmidtfTable& dmidtfs,
            const phaseSystem::dmdtfTable& Tfs,
            const scalar weight,
            const latentHeatScheme scheme,
            phaseSystem::heatTransferTable& eqns
        ) const;

        //- Add energy transfer terms which result from specie phase changes
        //  Weight is the proportion of the latent heat contribution applied to
        //  the downwind side of the mass transfer.
        void addDmidtHefs
        (
            const phaseSystem::dmidtfTable& dmidtfs,
            const phaseSystem::dmdtfTable& Tfs,
            const scalar weight,
            const latentHeatScheme scheme,
            phaseSystem::heatTransferTable& eqns
        ) const;


    // Protected data

        //- Residual mass fraction used for linearisation of heat transfer
        //  terms that result from mass transfers of individual species
        scalar residualY_;



public:

    // Constructors

        //- Construct from fvMesh
        HeatTransferPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~HeatTransferPhaseSystem();


    // Member Functions

        //- Return the latent heat for a given pair, mass transfer rate (used
        //  only for it's sign), and interface temperature
        virtual tmp<volScalarField> L
        (
            const phasePair& pair,
            const volScalarField& dmdtf,
            const volScalarField& Tf,
            const latentHeatScheme scheme
        ) const;

        //- As above, but for a cell-set
        virtual tmp<scalarField> L
        (
            const phasePair& pair,
            const scalarField& dmdtf,
            const scalarField& Tf,
            const labelUList& cells,
            const latentHeatScheme scheme
        ) const;

        //- Return the latent heat for a given pair, specie, mass transfer rate
        //  (used only for it's sign), and interface temperature
        virtual tmp<volScalarField> Li
        (
            const phasePair& pair,
            const word& member,
            const volScalarField& dmidtf,
            const volScalarField& Tf,
            const latentHeatScheme scheme
        ) const;

        //- As above, but for a cell-set
        virtual tmp<scalarField> Li
        (
            const phasePair& pair,
            const word& member,
            const scalarField& dmidtf,
            const scalarField& Tf,
            const labelUList& cells,
            const latentHeatScheme scheme
        ) const;

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HeatTransferPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
