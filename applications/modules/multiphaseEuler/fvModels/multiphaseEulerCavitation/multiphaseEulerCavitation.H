/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::fv::multiphaseEulerCavitation

Description
    Model for cavitation phase change between two phases.

    Example usage:
    \verbatim
    cavitation
    {
        type            multiphaseEuler::cavitation;
        libs            ("libmultiphaseEulerFvModels.so");

        // Note: Order is important. The liquid should be specified first.
        phases          (liquid gas);

        energySemiImplicit no;

        model           Kunz;

        pSat            79995.75943;

        UInf            5.33;
        tInf            0.028142589;
        Cc              100;
        Cv              100;
    }
    \endverbatim

SourceFiles
    multiphaseEulerCavitation.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseEulerCavitation_H
#define multiphaseEulerCavitation_H

#include "phaseChange.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
namespace compressible
{
    class cavitationModel;
}

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class multiphaseEulerCavitation Declaration
\*---------------------------------------------------------------------------*/

class multiphaseEulerCavitation
:
    public phaseChange
{
private:

    // Private Data

        //- Reference to the phase system
        const phaseSystem& fluid_;

        //- Reference to the liquid phase
        const phaseModel& liquid_;

        //- Reference to the vapour phase
        const phaseModel& vapour_;

        //- Interface
        const phaseInterface interface_;

        //- Reference to the field associated with the pressure equation
        const volScalarField& p_rgh_;

        //- Cavitation model
        autoPtr<compressible::cavitationModel> cavitationModel_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("multiphaseEuler::cavitation");


    // Constructors

        //- Construct from explicit source name and mesh
        multiphaseEulerCavitation
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // Evaluation

            //- Return the fraction of the latent heat that is transferred into
            //  the second phase
            virtual tmp<DimensionedField<scalar, volMesh>> Lfraction() const;


        // Sources

            //- Return the mass transfer rate
            virtual tmp<DimensionedField<scalar, volMesh>> mDot() const;

            //- Use phaseChange's source functions
            using phaseChange::addSup;

            //- Override the pressure equation to add the mass transfer rate
            //  linearised in the pressure
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn
            ) const;


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
