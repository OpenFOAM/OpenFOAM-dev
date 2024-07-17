/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::fv::constantIgnition

Description
    Simple constant rate ignition model for the Weller b-Xi combustion models

Usage
    Example usage:
    \verbatim
    constantIgnition1
    {
        type            constantIgnition;

        select          cellZone;
        cellZone        ignition;

        duration        0.003;
        strength        2;
    }
    \endverbatim

    Where:
    \table
        Property     | Description             | Required    | Default value
        duration     | Ignition duration       | yes         |
        strength     | Ignition strength [1/s] | yes         |
    \endtable

SourceFiles
    constantIgnition.C

\*---------------------------------------------------------------------------*/

#ifndef constantIgnition_H
#define constantIgnition_H

#include "fvModel.H"
#include "fvCellSet.H"
#include "XiFluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class constantIgnition Declaration
\*---------------------------------------------------------------------------*/

class constantIgnition
:
    public fvModel
{
    // Private Data

        //- The XiFluid solver
        const solvers::XiFluid& XiFluid_;

        //- The set of ignition cells
        fvCellSet set_;

        //- Ignition duration
        scalar duration_;

        //- Ignition strength
        scalar strength_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();


public:

    //- Runtime type information
    TypeName("constantIgnition");


    // Constructors

        //- Construct from explicit source name and mesh
        constantIgnition
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        constantIgnition
        (
            const constantIgnition&
        ) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Add explicit and implicit contributions to compressible equation

            //- Add ignition contribution to b equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& b,
                fvMatrix<scalar>& eqn
            ) const;


        // Mesh motion

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const constantIgnition&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
