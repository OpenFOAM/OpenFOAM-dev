/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::fv::interfaceTurbulenceDamping

Description
    Free-surface phase turbulence damping function

    Adds an extra source term to the mixture or phase epsilon or omega
    equation to reduce turbulence generated near a free-surface.  The
    implementation is based on

    Reference:
    \verbatim
        Frederix, E. M. A., Mathur, A., Dovizio, D., Geurts, B. J.,
        & Komen, E. M. J. (2018).
        Reynolds-averaged modeling of turbulence damping
        near a large-scale interface in two-phase flow.
        Nuclear engineering and design, 333, 122-130.
    \endverbatim

    but with an improved formulation for the coefficient \c A appropriate for
    unstructured meshes including those with split-cell refinement patterns.
    However the dimensioned length-scale coefficient \c delta remains and must
    be set appropriately for the case by performing test runs and comparing with
    known results.  Clearly this model is far from general and more research is
    needed in order that \c delta can be obtained directly from the interface
    flow and turbulence conditions.

Usage
    Example usage:
    \verbatim
    interfaceTurbulenceDamping
    {
        type    interfaceTurbulenceDamping;

        libs    ("libmultiphaseEulerFoamFvModels.so");

        phase   water;

        // Interface turbulence damping length scale
        // This is a required input as described in section 3.3 of the paper
        delta   1e-4;
    }
    \endverbatim

SourceFiles
    interfaceTurbulenceDamping.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceTurbulenceDamping_H
#define interfaceTurbulenceDamping_H

#include "fvModel.H"
#include "phaseSystem.H"
#include "phaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class interfaceTurbulenceDamping Declaration
\*---------------------------------------------------------------------------*/

class interfaceTurbulenceDamping
:
    public fvModel
{
    // Private Data

        //- The name of the Lagrangian phase
        word phaseName_;

        //- Field name
        word fieldName_;

        //- Interface turbulence damping length scale
        //  This is a required input as described in section 3.3 of the paper
        dimensionedScalar delta_;

        //- Reference to the phase
        const phaseModel& phase_;

        //- Reference to the mixture turbulence model
        const phaseCompressible::momentumTransportModel& turbulence_;

        // Turbulence model coefficients

            dimensionedScalar C2_;
            dimensionedScalar betaStar_;
            dimensionedScalar beta_;


    // Private Member Functions

        //- Interface fraction in a cell
        tmp<volScalarField::Internal> interfaceFraction
        (
            const volScalarField& alpha
        ) const;

        //- Add explicit contribution to incompressible or compressible
        //  mixture epsilon or omega equation
        template<class RhoType>
        void addRhoSup
        (
            const RhoType& rho,
            fvMatrix<scalar>& eqn,
            const word& fieldName
        ) const;


public:

    //- Runtime type information
    TypeName("interfaceTurbulenceDamping");


    // Constructors

        //- Construct from explicit source name and mesh
        interfaceTurbulenceDamping
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        interfaceTurbulenceDamping
        (
            const interfaceTurbulenceDamping&
        ) = delete;


    // Member Functions

        //- Return the list of fields for which the option adds source term
        //  to the transport equation
        virtual wordList addSupFields() const;

        //- Add explicit contribution to mixture epsilon or omega equation
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const word& fieldName
        ) const;

        //- Add explicit contribution to compressible
        //  mixture epsilon or omega equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const word& fieldName
        ) const;

        //- Add explicit contribution to phase epsilon or omega equation
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const word& fieldName
        ) const;


        // Mesh changes

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
        void operator=(const interfaceTurbulenceDamping&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
