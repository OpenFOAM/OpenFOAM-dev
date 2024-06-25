/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2024 OpenFOAM Foundation
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
    Foam::fv::buoyancyEnergy

Description
    Calculates and applies the buoyancy energy source rho*(U&g) to the energy
    equation.

Usage
    Example usage:
    \verbatim
    buoyancyEnergy1
    {
        type        buoyancyEnergy;

        U           U;  // Name of the velocity field
    }
    \endverbatim

SourceFiles
    buoyancyEnergy.C

\*---------------------------------------------------------------------------*/

#ifndef buoyancyEnergy_H
#define buoyancyEnergy_H

#include "fvModel.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class buoyancyEnergy Declaration
\*---------------------------------------------------------------------------*/

class buoyancyEnergy
:
    public fvModel
{
    // Private Data

        //- Optional phase name
        word phaseName_;

        //- Name of velocity field; default = U.<phase>
        word UName_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();


public:

    //- Runtime type information
    TypeName("buoyancyEnergy");


    // Constructors

        //- Construct from explicit source name and mesh
        buoyancyEnergy
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        buoyancyEnergy(const buoyancyEnergy&) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Evaluate

            //- Add explicit contribution to compressible energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add explicit contribution to phase energy equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const buoyancyEnergy&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
