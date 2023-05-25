/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    Foam::fv::VoFClouds

Description
    Lagrangian clouds model for VoF simulations.

Usage
    Example usage:
    \verbatim
    VoFClouds
    {
        type    VoFClouds;
        phase   water;
    }
    \endverbatim

SourceFiles
    VoFClouds.C

\*---------------------------------------------------------------------------*/

#ifndef VoFClouds_H
#define VoFClouds_H

#include "fvModel.H"
#include "fluidThermo.H"
#include "parcelCloudList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class VoFClouds Declaration
\*---------------------------------------------------------------------------*/

class VoFClouds
:
    public fvModel
{
    // Private Data

        //- The name of the Lagrangian phase
        word phaseName_;

        //- The name of the Lagrangian carrier phase
        word carrierPhaseName_;

        //- Reference to the Lagrangian thermo
        const fluidThermo& thermo_;

        //- Reference to the carrier phase thermo
        const fluidThermo& carrierThermo_;

        //- The Lagrangian cloud list
        mutable parcelCloudList clouds_;

        //- Current time index (used for updating)
        mutable label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("VoFClouds");


    // Constructors

        //- Construct from explicit source name and mesh
        VoFClouds
        (
            const word& sourceName,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        VoFClouds
        (
            const VoFClouds&
        ) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Correct

            //- Solve the Lagrangian clouds and update the sources
            virtual void correct();


        // Add explicit and implicit contributions to compressible equation

            //- Add explicit contribution to phase continuity
            virtual void addSup
            (
                const volScalarField& alpha,
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

            //- Add explicit contribution to phase energy equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

            //- Add implicit contribution to mixture momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;


        // Mesh changes

            //- Prepare for mesh_.update
            virtual void preUpdateMesh();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const VoFClouds&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
