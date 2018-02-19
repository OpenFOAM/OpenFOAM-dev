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
    Foam::HeatAndMassTransferPhaseSystem

Description
    Base class to support interfacial heat and mass transfer between a number
    of phases. Can be used on its own to simulate flows without mass transfer
    with two-resistance interfacial heat trasnfer.

SourceFiles
    HeatAndMassTransferPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef HeatAndMassTransferPhaseSystem_H
#define HeatAndMassTransferPhaseSystem_H

#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class modelType>
class BlendedInterfacialModel;

class blendingMethod;
class heatTransferModel;
class massTransferModel;

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
                autoPtr<BlendedInterfacialModel<heatTransferModel>>
            >,
            phasePairKey,
            phasePairKey::hash
        > heatTransferModelTable;

        typedef HashTable
        <
            HashTable
            <
                autoPtr<BlendedInterfacialModel<massTransferModel>>
            >,
            phasePairKey,
            phasePairKey::hash
        > massTransferModelTable;


    // Protected data

        //- Interface temperatures
        HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash> Tf_;

        // Sub Models

            //- Heat transfer models
            heatTransferModelTable heatTransferModels_;

            //- Mass transfer models
            massTransferModelTable massTransferModels_;


public:

    // Constructors

        //- Construct from fvMesh
        HeatAndMassTransferPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~HeatAndMassTransferPhaseSystem();


    // Member Functions

        //- Return true if there is mass transfer
        virtual bool transfersMass() const;

        //- Return the interfacial mass flow rate
        virtual tmp<volScalarField> dmdt(const phasePairKey& key) const;

        //- Return the total interfacial mass transfer rate for phase
        virtual tmp<volScalarField> dmdt(const phaseModel& phase) const;

        //- Return the momentum transfer matrices
        virtual autoPtr<phaseSystem::momentumTransferTable>
            momentumTransfer() const;

        //- Return the momentum transfer matrices for the face based algorithm
        virtual autoPtr<phaseSystem::momentumTransferTable>
            momentumTransferf() const;

        //- Return the heat transfer matrices
        virtual autoPtr<phaseSystem::heatTransferTable>
            heatTransfer() const;

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
    #include "HeatAndMassTransferPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
