/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::fv::wallPhaseChange

Description
    Base class for fvModels that represent phase change at a wall

See also
    Foam::alphatPhaseChangeWallFunctionFvPatchScalarField

SourceFiles
    wallPhaseChange.C

\*---------------------------------------------------------------------------*/

#ifndef wallPhaseChange_H
#define wallPhaseChange_H

#include "phaseChange.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class wallPhaseChange Declaration
\*---------------------------------------------------------------------------*/

class wallPhaseChange
:
    public phaseChange
{
    // Private Data

        //- Reference to the phase system
        const phaseSystem& fluid_;

        //- References to the phases
        const Pair<const phaseModel&> phases_;

        //- References to the turbulent thermal diffusivities
        const Pair<const volScalarField&> alphats_;


protected:

    // Protected Member Functions

        //- Return boundary conditions for the phase change rate field. Use the
        //  specified active type for all phase-changing boundaries, and a
        //  zero-fixed-value type for others.
        wordList mDotBoundaryTypes(const word& activeType) const;


public:

    //- Runtime type information
    TypeName("wallPhaseChange");


    // Constructors

        //- Construct from explicit source name and mesh
        wallPhaseChange
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const wordList& species
        );


    // Member Functions

        // Access

            //- Access the phase system
            inline const phaseSystem& fluid() const;

            //- Access the phases
            inline const Pair<const phaseModel&>& phases() const;

            //- Access the turbulent thermal diffusivities
            inline const Pair<const volScalarField&>& alphats() const;

            //- Return whether phase change can occur on the given patch
            bool isPatchActive(const label patchi) const;

            //- Return a mask indicating whether phase change is occurring
            virtual const scalarField& active(const label patchi) const = 0;

            //- Access the turbulent thermal diffusivities for a patch
            virtual Pair<const scalarField&> alphats(const label) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallPhaseChangeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
