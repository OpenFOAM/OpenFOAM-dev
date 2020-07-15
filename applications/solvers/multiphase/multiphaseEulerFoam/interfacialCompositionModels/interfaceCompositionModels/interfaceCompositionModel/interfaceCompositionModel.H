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
    Foam::interfaceCompositionModel

Description
    Generic base class for interface composition models. These models describe
    the composition in phase 1 of the supplied pair at the interface with phase
    2.

SourceFiles
    interfaceCompositionModel.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceCompositionModel_H
#define interfaceCompositionModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "hashedWordList.H"
#include "rhoReactionThermo.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phaseModel;
class phasePair;

/*---------------------------------------------------------------------------*\
                        Class interfaceCompositionModel Declaration
\*---------------------------------------------------------------------------*/

class interfaceCompositionModel
{
    // Private data

        //- Phase pair
        const phasePair& pair_;

        //- Names of the transferring species
        const hashedWordList species_;

        //- Lewis number
        const dimensionedScalar Le_;

        //- Phase
        const phaseModel& phase_;

        //- Other phase
        const phaseModel& otherPhase_;

        //- Multi-component thermo model for this side of the interface
        const rhoReactionThermo& thermo_;

        //- General thermo model for the other side of the interface
        const rhoThermo& otherThermo_;


public:

    //- Runtime type information
    TypeName("interfaceCompositionModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            interfaceCompositionModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair
            ),
            (dict, pair)
        );


    // Constructors

        //- Construct from a dictionary and a phase pair
        interfaceCompositionModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~interfaceCompositionModel();


    // Selectors

        static autoPtr<interfaceCompositionModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        // Access

            //- Return the phase pair
            inline const phasePair& pair() const;

            //- Return the transferring species names
            inline const hashedWordList& species() const;

            //- Return the thermo
            inline const rhoReactionThermo& thermo() const;

            //- Return the composition
            inline const basicSpecieMixture& composition() const;

            //- Return the other thermo
            inline const rhoThermo& otherThermo() const;

            //- Return whether the other side has a multi-specie composition
            inline bool otherHasComposition() const;

            //- Return the other composition
            inline const basicSpecieMixture& otherComposition() const;


        // Evaluation

            //- Interface mass fraction
            virtual tmp<volScalarField> Yf
            (
                const word& speciesName,
                const volScalarField& Tf
            ) const = 0;

            //- The interface mass fraction derivative w.r.t. temperature
            virtual tmp<volScalarField> YfPrime
            (
                const word& speciesName,
                const volScalarField& Tf
            ) const = 0;

            //- Mass fraction difference between the interface and the field
            virtual tmp<volScalarField> dY
            (
                const word& speciesName,
                const volScalarField& Tf
            ) const;

            //- Mass diffusivity
            virtual tmp<volScalarField> D
            (
                const word& speciesName
            ) const;

            //- Latent heat
            virtual tmp<volScalarField> L
            (
                const word& speciesName,
                const volScalarField& Tf
            ) const;

            //- Add latent heat flow rate to total
            virtual void addDmdtL
            (
                const volScalarField& K,
                const volScalarField& Tf,
                volScalarField& dmdtL,
                volScalarField& dmdtLPrime
            ) const;


        //- Update the composition
        virtual void update(const volScalarField& Tf) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interfaceCompositionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
