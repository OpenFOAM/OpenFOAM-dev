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
    Foam::fv::singleComponentPhaseChangeBase

Description
    Base class for phase change models in which only a single component changes
    phase. Can be applied to any combination of pure and multicomponent phases.
    If either phase is multicomponent, then a single specie must be identified
    as the one that changes phase.

SourceFiles
    singleComponentPhaseChangeBase.C

\*---------------------------------------------------------------------------*/

#ifndef singleComponentPhaseChangeBase_H
#define singleComponentPhaseChangeBase_H

#include "phaseChangeBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class singleComponentPhaseChangeBase Declaration
\*---------------------------------------------------------------------------*/

class singleComponentPhaseChangeBase
:
    public phaseChangeBase
{
private:

    // Private Data

        //- The name of the changing specie, or word::null if neither thermo
        //  is multicomponent
        const word specie_;

        //- The indices of the changing species in the phases, or -1 if the
        //  phases are not multicomponent
        const Pair<label> specieis_;

        //- Whether or not to linearise the energy source
        bool energySemiImplicit_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();


public:

    //- Runtime type information
    TypeName("singleComponentPhaseChangeBase");


    // Constructors

        //- Construct from explicit source name and mesh
        singleComponentPhaseChangeBase
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const Pair<bool>& fluidThermosRequired,
            const Pair<bool>& specieThermosRequired
        );


    // Member Functions

        // Access

            //- Return the name of the changing specie
            inline const word& specie() const;

            //- Return the indices of the changing species in the phases
            inline const Pair<label>& specieis() const;


        // Sources

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

#include "singleComponentPhaseChangeBaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
