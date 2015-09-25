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
    Foam::ThermalPhaseChangePhaseSystem

Description
    Class to provide interfacial heat and mass transfer between a number of
    phases according the interfacial temperature approximated by the saturation
    temperature.

    Currently only a single specified specie is considered volatile and changes
    phase, all other species are considered nonvolatile and do not
    affect the mass-transfer.

SourceFiles
    ThermalPhaseChangePhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef ThermalPhaseChangePhaseSystem_H
#define ThermalPhaseChangePhaseSystem_H

#include "HeatAndMassTransferPhaseSystem.H"
#include "saturationModel.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class ThermalPhaseChangePhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class ThermalPhaseChangePhaseSystem
:
    public HeatAndMassTransferPhaseSystem<BasePhaseSystem>
{

protected:

    // Protected data

        //- Name of the volatile specie
        word volatile_;

        //- The saturation model used to evaluate Tsat = Tf
        autoPtr<saturationModel> saturationModel_;

        // Mass transfer enabled
        Switch massTransfer_;

        //- Interfacial Mass transfer rate
        HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>
            iDmdt_;

        //- Wall Mass transfer rate
        HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>
            wDmdt_;


public:

    // Constructors

        //- Construct from fvMesh
        ThermalPhaseChangePhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~ThermalPhaseChangePhaseSystem();


    // Member Functions

        //- Return the saturationModel
        const saturationModel& saturation() const;

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
#   include "ThermalPhaseChangePhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
