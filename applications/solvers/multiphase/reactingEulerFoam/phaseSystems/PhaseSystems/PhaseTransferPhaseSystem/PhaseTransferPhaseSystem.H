/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::MassTransferTransferPhaseSystem

Description
    Class which models non-thermally-coupled mass transfers; i.e.,
    representation changes, rather than phase changes.

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
protected:

    // Protected typedefs

        typedef HashTable
        <
            autoPtr<BlendedInterfacialModel<phaseTransferModel>>,
            phasePairKey,
            phasePairKey::hash
        > phaseTransferModelTable;

        typedef HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>
            rDmdtTable;


    // Protected data

        // Sub Models

            //- Mass transfer models
            phaseTransferModelTable phaseTransferModels_;

            //- Mass transfer rates
            rDmdtTable rDmdt_;

    // Protected member functions

        //- Return the representation mass transfer rate
        virtual tmp<volScalarField> rDmdt(const phasePairKey& key) const;


public:

    // Constructors

        //- Construct from fvMesh
        PhaseTransferPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~PhaseTransferPhaseSystem();


    // Member Functions

        //- Return the mass transfer rate for a pair
        virtual tmp<volScalarField> dmdt(const phasePairKey& key) const;

        //- Return the mass transfer rates for each phase
        virtual Xfer<PtrList<volScalarField>> dmdts() const;

        //- Return the mass transfer matrices
        virtual autoPtr<phaseSystem::massTransferTable> massTransfer() const;

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
