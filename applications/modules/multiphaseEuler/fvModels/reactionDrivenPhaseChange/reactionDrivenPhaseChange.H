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
    Foam::fv::reactionDrivenPhaseChange

Description
    Model for mass-diffusion rate limited phase change between two phases.

    This model will transfer species out of a reacting phase at exactly the
    rate at which the reactions produce them. It facilitates a rough
    approximation to surface reactions for which certain species get consumed
    at the interface.

Usage
    Example usage:
    \verbatim
    phaseChange
    {
        type            reactionDrivenPhaseChange;
        libs            ("libmultiphaseEulerFvModels.so");

        phases          (particles vapour);

        species         (TiO2);

        energySemiImplicit no;
    }
    \endverbatim

SourceFiles
    reactionDrivenPhaseChange.C

\*---------------------------------------------------------------------------*/

#ifndef reactionDrivenPhaseChange_H
#define reactionDrivenPhaseChange_H

#include "phaseChange.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
              Class reactionDrivenPhaseChange Declaration
\*---------------------------------------------------------------------------*/

class reactionDrivenPhaseChange
:
    public phaseChange
{
private:

    // Private Data

        //- Phase system
        const phaseSystem& fluid_;

        //- Phase 1
        const phaseModel& phase1_;

        //- Phase 2
        const phaseModel& phase2_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("reactionDrivenPhaseChange");


    // Constructors

        //- Construct from explicit source name and mesh
        reactionDrivenPhaseChange
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // Sources

            //- Return the mass transfer rate of a specie
            virtual tmp<DimensionedField<scalar, volMesh>> mDot
            (
                const label mDoti
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

// ************************************************************************* //
