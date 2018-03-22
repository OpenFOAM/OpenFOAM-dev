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
    Foam::HeatTransferPhaseSystem

Description
    Class which models interfacial heat transfer between a number of phases.

SourceFiles
    HeatTransferPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef HeatTransferPhaseSystem_H
#define HeatTransferPhaseSystem_H

#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blendingMethod;
template<class modelType> class BlendedInterfacialModel;
class heatTransferModel;

/*---------------------------------------------------------------------------*\
                    Class HeatTransferPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class HeatTransferPhaseSystem
:
    public BasePhaseSystem
{
protected:

    // Protected typedefs

        typedef HashTable
        <
            autoPtr<BlendedInterfacialModel<heatTransferModel>>,
            phasePairKey,
            phasePairKey::hash
        > heatTransferModelTable;


    // Protected data

        // Sub Models

            //- Heat transfer models
            heatTransferModelTable heatTransferModels_;


public:

    // Constructors

        //- Construct from fvMesh
        HeatTransferPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~HeatTransferPhaseSystem();


    // Member Functions

        //- Return the heat transfer matrices
        virtual autoPtr<phaseSystem::heatTransferTable> heatTransfer() const;

        //- Return the mass transfer matrices
        virtual autoPtr<phaseSystem::massTransferTable> massTransfer() const;

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
