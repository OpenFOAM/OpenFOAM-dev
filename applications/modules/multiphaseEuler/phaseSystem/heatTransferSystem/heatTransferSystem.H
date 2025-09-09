/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2025 OpenFOAM Foundation
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
    Foam::heatTransferSystem

Description
    Class which provides interfacial momentum transfer between a number of
    phases. Drag, virtual mass, lift, wall lubrication and turbulent dispersion
    are all modelled. Access is provided to a various forces and coefficients
    needed by the solver.

SourceFiles
    heatTransferSystem.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferSystem_H
#define heatTransferSystem_H

#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class blendedHeatTransferModel;
class blendedSidedHeatTransferModel;

/*---------------------------------------------------------------------------*\
                 Class heatTransferSystem Declaration
\*---------------------------------------------------------------------------*/

class heatTransferSystem
:
    public IOdictionary
{
    // Private Type Definitions

        //- Table of one-resistance heat transfer models
        typedef
            HashPtrTable
            <
                blendedHeatTransferModel,
                phaseInterfaceKey,
                phaseInterfaceKey::hash
            > modelsTable;

        //- Table of two-resistance heat transfer models
        typedef
            HashPtrTable
            <
                blendedSidedHeatTransferModel,
                phaseInterfaceKey,
                phaseInterfaceKey::hash
            > sidedModelsTable;


    // Private Data

        //- Reference to the phase system
        const phaseSystem& fluid_;

        //- One-resistance heat transfer models
        modelsTable models_;

        //- Two-resistance heat transfer models
        sidedModelsTable sidedModels_;


    // Private Member Functions

        //- Create IO object for an optional constant/heatTransfer file
        IOobject io(const phaseSystem&) const;

        //- Return the models dictionary
        const dictionary& modelsDict() const;

        //- Return the two-resistance heat transfer coefficients between a pair
        //  of phases with additional arguments
        template<class ... Args>
        Pair<tmp<volScalarField>> Hs
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            Args ... args
        ) const;


public:

    //- Runtime type information
    TypeName("heatTransferSystem");


    //- Default name of the phase properties dictionary
    static const word propertiesName;


    // Constructors

        //- Construct from a phase system
        heatTransferSystem(const phaseSystem&);


    //- Destructor
    virtual ~heatTransferSystem();


    // Member Functions

        //- Return the two-resistance heat transfer coefficients between a pair
        //  of phases
        Pair<tmp<volScalarField>> Hs
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Return the two-resistance heat transfer coefficients between a pair
        //  of phases with a specified residual volume fraction
        Pair<tmp<volScalarField>> Hs
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const scalar residualAlpha
        ) const;

        //- Return the heat transfer matrices
        autoPtr<HashPtrTable<fvScalarMatrix>> heatTransfer() const;

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
