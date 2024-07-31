/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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
    Foam::fv::multicomponentPhaseChange

Description
    Base class for phase change models in which multiple components may change
    phase. This can only be applied between multicomponent phases.

SourceFiles
    multicomponentPhaseChange.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentPhaseChange_H
#define multicomponentPhaseChange_H

#include "phaseChange.H"
#include "hashedWordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class multicomponentPhaseChange Declaration
\*---------------------------------------------------------------------------*/

class multicomponentPhaseChange
:
    public phaseChange
{
private:

    // Private Data

        //- The names of the changing species
        const hashedWordList species_;

        //- Whether or not to linearise the energy source
        bool energySemiImplicit_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("multicomponentPhaseChange");


    // Constructors

        //- Construct from explicit source name and mesh
        multicomponentPhaseChange
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const Pair<bool>& fluidThermosRequired
        );


    // Member Functions

        // Access

            //- Return the names of the changing species
            inline const hashedWordList& species() const;


        // Sources

            //- Return the mass transfer rates of each specie
            virtual tmp<DimensionedField<scalar, volMesh>> mDot
            (
                const label mDoti
            ) const = 0;

            //- Return the total mass transfer rate
            virtual tmp<DimensionedField<scalar, volMesh>> mDot() const;

            //- Override the energy equation to add the phase change heat, or
            //  the species equations to add the relevant mass sources
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& heOrYi,
                fvMatrix<scalar>& eqn
            ) const;


        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentPhaseChangeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
