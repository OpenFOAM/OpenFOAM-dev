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
    Foam::fv::clouds

Description
    This fvModel adds any number of Lagrangian clouds to any single-phase
    solver. The particles are tracked through, and exchange sources with, the
    Eulerian flow field.

    As well as the fvModel controls, properties must be specified for each
    cloud. For a single cloud, these should be provided in the
    constant/cloudProperties file. For multiple clouds, the list of cloud names
    must first be provided in the constant/clouds file. Then, each named cloud
    has its properties specified in its constant/\<cloudName\>Properties file.

    The application of sources to the Eulerian fields is controlled by the
    solution/coupled switch in each cloud's properties file. If set to "true"
    then the Eulerian phase will have forces, and heat and mass sources applied
    to it by the Lagrangian phase. If set to "false" then these will be omitted,
    and the Lagrangian phase will not affect the Eulerian phase.

    If this model is used with an incompressible solver, then the density of
    the Eulerian phase must be specified in the constant/physicalProperties
    dictionary.

    Gravity will be read from the constant/g file if present, otherwise it will
    default to zero.

Usage
    Example usage:
    \verbatim
    clouds
    {
        libs        ("liblagrangianParcel.so");
        type        clouds;
    }
    \endverbatim

    \table
        Property | Description                     | Required   | Default value
        type     | Type name: clouds               | yes        |
        rho      | Name of the density field       | no         | rho
        U        | Name of the velocity field      | no         | U
    \endtable

SourceFiles
    clouds.C

\*---------------------------------------------------------------------------*/

#ifndef clouds_H
#define clouds_H

#include "fvModel.H"
#include "fluidThermo.H"
#include "viscosityModel.H"
#include "uniformDimensionedFields.H"
#include "parcelCloudList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class clouds Declaration
\*---------------------------------------------------------------------------*/

class clouds
:
    public fvModel
{
    // Private Data

        //- Optional acceleration due to gravity
        const uniformDimensionedVectorField g_;

        //- Flag to indicate whether a thermo model is available for the
        //  carrier
        const bool carrierHasThermo_;

        //- Reference to the carrier phase thermo. Valid only if the carrier
        //  has thermo.
        const tmpNrc<fluidThermo> tCarrierThermo_;

        //- Reference to the carrier viscosity model. Valid only if the carrier
        //  does not have thermo.
        const tmpNrc<viscosityModel> tCarrierViscosity_;

        //- Density field. Valid only if the carrier does not have thermo.
        const tmp<volScalarField> tRho_;

        //- Viscosity field. Valid only if the carrier does not have thermo.
        tmp<volScalarField> tMu_;

        //- Names of the clouds
        const wordList cloudNames_;

        //- Name of the density field
        const word rhoName_;

        //- Name of the velocity field
        const word UName_;

        //- The Lagrangian cloud list
        mutable autoPtr<parcelCloudList> cloudsPtr_;

        //- Current time index (used for updating)
        mutable label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("clouds");


    // Constructors

        //- Construct from explicit source name and mesh
        clouds
        (
            const word& sourceName,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        clouds
        (
            const clouds&
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

            //- Add source to continuity equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add source to enthalpy or species equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& field,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add source to incompressible momentum equation
            virtual void addSup
            (
                const volVectorField& U,
                fvMatrix<vector>& eqn
            ) const;

            //- Add source to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volVectorField& U,
                fvMatrix<vector>& eqn
            ) const;


        // Mesh changes

            //- Prepare for mesh update
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
        void operator=(const clouds&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
