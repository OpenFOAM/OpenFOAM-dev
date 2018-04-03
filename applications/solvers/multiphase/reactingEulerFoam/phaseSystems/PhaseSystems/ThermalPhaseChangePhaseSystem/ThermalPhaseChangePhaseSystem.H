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
    Foam::ThermalPhaseChangePhaseSystem

Description
    Class to provide interfacial heat and mass transfer between a number of
    phases according the interfacial temperature approximated by the saturation
    temperature.

SourceFiles
    ThermalPhaseChangePhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef ThermalPhaseChangePhaseSystem_H
#define ThermalPhaseChangePhaseSystem_H

#include "phaseSystem.H"
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
    public BasePhaseSystem
{

protected:

    // Protected typedefs

        typedef HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>
            iDmdtTable;

        typedef HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>
            wDmdtTable;

        typedef HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>
            wMDotLTable;


    // Protected data

        //- Name of the volatile specie
        word volatile_;

        //- The saturation model used to evaluate Tsat = Tf
        autoPtr<saturationModel> saturationModel_;

        // Phase change enabled
        Switch phaseChange_;

        //- Interfacial Mass transfer rate
        iDmdtTable iDmdt_;

        //- Boundary Mass transfer rate
        wDmdtTable wDmdt_;

        //- Boundary thermal energy transfer rate
        wMDotLTable wMDotL_;


    // Protected member functions

        //- Return the interfacial mass transfer rate for a pair
        tmp<volScalarField> iDmdt(const phasePairKey& key) const;

        //- Return the boundary mass transfer rate for a pair
        tmp<volScalarField> wDmdt(const phasePairKey& key) const;


public:

    // Constructors

        //- Construct from fvMesh
        ThermalPhaseChangePhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~ThermalPhaseChangePhaseSystem();


    // Member Functions

        //- Return the saturationModel
        const saturationModel& saturation() const;

        //- Return the mass transfer rate for a pair
        virtual tmp<volScalarField> dmdt(const phasePairKey& key) const;

        //- Return the mass transfer rates for each phase
        virtual Xfer<PtrList<volScalarField>> dmdts() const;

        //- Return the heat transfer matrices
        virtual autoPtr<phaseSystem::heatTransferTable> heatTransfer() const;

        //- Return the mass transfer matrices
        virtual autoPtr<phaseSystem::massTransferTable> massTransfer() const;

        //- Correct the interface thermodynamics
        virtual void correctInterfaceThermo();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ThermalPhaseChangePhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
