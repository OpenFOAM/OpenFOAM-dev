/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::HeatAndMassTransferPhaseSystem

Description
    Class which models interfacial heat and mass transfer between a number of
    phases. Mass is transferred to or from a phase according to a composition
    model at the surface of another phase. Heat is transferred from a phase to
    an interfacial temperature. The interface temperature is calculated such
    that the net rate at which the heat is transferred to the interface is
    equal to the latent heat consumed by the mass transfer.

SourceFiles
    HeatAndMassTransferPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef HeatAndMassTransferPhaseSystem_H
#define HeatAndMassTransferPhaseSystem_H

#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template <class modelType>
class BlendedInterfacialModel;

class blendingMethod;
class heatTransferModel;
class massTransferModel;
class interfaceCompositionModel;

/*---------------------------------------------------------------------------*\
                 Class HeatAndMassTransferPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class HeatAndMassTransferPhaseSystem
:
    public BasePhaseSystem
{
protected:

    // Protected typedefs

        typedef HashTable
        <
            HashTable
            <
                autoPtr<BlendedInterfacialModel<heatTransferModel> >
            >,
            phasePairKey,
            phasePairKey::hash
        > heatTransferModelTable;

        typedef HashTable
        <
            HashTable
            <
                autoPtr<BlendedInterfacialModel<massTransferModel> >
            >,
            phasePairKey,
            phasePairKey::hash
        > massTransferModelTable;

        typedef HashTable
        <
            autoPtr<interfaceCompositionModel>,
            phasePairKey,
            phasePairKey::hash
        > interfaceCompositionModelTable;


    // Protected data

        //- Mass transfer rate
        HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>
            dmdt_;

        //- Explicit part of the mass transfer rate
        HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>
            dmdtExplicit_;

        //- Interface temperatures
        HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash> Tf_;

        // Sub Models

            //- Heat transfer models
            heatTransferModelTable heatTransferModels_;

            //- Mass transfer models
            massTransferModelTable massTransferModels_;

            //- Interface composition models
            interfaceCompositionModelTable interfaceCompositionModels_;


public:

    // Constructors

        //- Construct from fvMesh
        HeatAndMassTransferPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~HeatAndMassTransferPhaseSystem();


    // Member Functions

        //- Return the interfacial mass flow rate
        virtual tmp<volScalarField> dmdt
        (
            const phasePairKey& key
        ) const;

        //- Return the momentum transfer matrices
        virtual autoPtr<phaseSystem::momentumTransferTable>
            momentumTransfer() const;

        //- Return the heat transfer matrices
        virtual autoPtr<phaseSystem::heatTransferTable>
            heatTransfer() const;

        //- Return the mass transfer matrices
        virtual autoPtr<phaseSystem::massTransferTable>
            massTransfer() const;

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "HeatAndMassTransferPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
