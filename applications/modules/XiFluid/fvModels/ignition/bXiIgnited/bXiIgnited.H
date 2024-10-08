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
    Foam::fv::bXiIgnited

Description
    Simple ignited model for the Weller b-Xi combustion models

    which sets ignited true and ignition false so that the combustion model
    executes continuously.

Usage
    Example usage:
    \verbatim
    bXiIgnited
    {
        type            bXiIgnited;
    }
    \endverbatim

SourceFiles
    bXiIgnited.C

\*---------------------------------------------------------------------------*/

#ifndef bXiIgnited_H
#define bXiIgnited_H

#include "bXiIgnition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class bXiIgnited Declaration
\*---------------------------------------------------------------------------*/

class bXiIgnited
:
    public bXiIgnition
{

public:

    //- Runtime type information
    TypeName("bXiIgnited");


    // Constructors

        //- Construct from explicit source name and mesh
        bXiIgnited
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        bXiIgnited
        (
            const bXiIgnited&
        ) = delete;


    // Member Functions

        // Checks

            //- Return true during the ignition duration
            virtual bool igniting() const;

            //- Return true during the combustion duration
            virtual bool ignited() const;


        // Add explicit and implicit contributions to compressible equation

            //- Add ignition contribution to b equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& b,
                fvMatrix<scalar>& eqn
            ) const;

            virtual void XiCorr
            (
                volScalarField& Xi,
                const volScalarField& b,
                const volScalarField& mgb
            ) const;


        // Mesh motion

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&)
            {}

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&)
            {}

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&)
            {}

            //- Update for mesh motion
            virtual bool movePoints()
            {
                return true;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const bXiIgnited&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
