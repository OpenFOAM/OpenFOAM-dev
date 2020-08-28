/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 OpenFOAM Foundation
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
    Foam::PhaseTransferPhaseSystem

Description
   Class which models non-thermally-coupled or weakly thermally coupled
   mass transfers.

SourceFiles
    PhaseTransferPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef PhaseTransferPhaseSystem_H
#define PhaseTransferPhaseSystem_H

#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blendingMethod;
template<class modelType> class BlendedInterfacialModel;
class phaseTransferModel;

/*---------------------------------------------------------------------------*\
                  Class PhaseTransferPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class PhaseTransferPhaseSystem
:
    public BasePhaseSystem
{
private:

    // Private typedefs

        typedef HashTable
        <
            autoPtr<BlendedInterfacialModel<phaseTransferModel>>,
            phasePairKey,
            phasePairKey::hash
        > phaseTransferModelTable;


    // Private data

        // Sub Models

            //- Mass transfer models
            phaseTransferModelTable phaseTransferModels_;

            //- Bulk mass transfer rates
            phaseSystem::dmdtfTable dmdtfs_;

            //- Specie mass transfer rates
            phaseSystem::dmidtfTable dmidtfs_;


    // Private Member Functions

        //- Return the summed mass transfers across each interface
        autoPtr<phaseSystem::dmdtfTable> totalDmdtfs() const;


protected:

    // Protected member functions

        //- Add specie transfer terms which result from bulk mass transfers
        void addDmdtYfs
        (
            const phaseSystem::dmdtfTable& dmdtfs,
            phaseSystem::specieTransferTable& eqns
        ) const;

        //- Add specie transfer terms which result from specie mass transfers
        void addDmidtYf
        (
            const phaseSystem::dmidtfTable& dmidtfs,
            phaseSystem::specieTransferTable& eqns
        ) const;


public:

    // Constructors

        //- Construct from fvMesh
        PhaseTransferPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~PhaseTransferPhaseSystem();


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

        //- Correct the mass transfer rates
        virtual void correct();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PhaseTransferPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
