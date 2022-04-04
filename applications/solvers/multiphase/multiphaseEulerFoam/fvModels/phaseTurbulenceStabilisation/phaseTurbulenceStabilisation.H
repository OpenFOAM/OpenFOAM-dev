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
    Foam::fv::phaseTurbulenceStabilisation

Description
    Phase turbulence stabilisation

    In the limit of a phase-fraction->0 the turbulence properties cannot be
    obtained from the phase turbulence model, coupling to the other phase/phases
    is required.  The phaseTurbulenceStabilisation fvModel stabilises the phase
    turbulence properties by adding transfer terms from the corresponding
    properties of the other phases when the phase-fraction is less than the
    specified \c alphaInversion.  This implementation is a generalisation of
    the approach used in the Foam::RASModels::LaheyKEpsilon and
    Foam::RASModels::continuousGasKEpsilon models to handle phase-inversion and
    free-surface flow and can be used with any combination of RAS turbulence
    models.

    To stabilise the solution of the phase turbulence equations \c
    alphaInversion can be set to a small value e.g. 1e-2, but unless the phase
    turbulence model is specifically designed to handle phase-inversion and both
    continuous and dispersed regimes it may be useful to set \c alphaInversion
    to a higher value, corresponding to the phase-fraction at which transision
    from continuous to dispersed happens and effectively use the turbulence
    properties of the other phase when the phase is dispersed.  This is of
    course an approximation to the real system and if accurate handling of both
    the continuous and dispersed phase regimes is required specially developed
    models should be used.

Usage
    Example usage:
    \verbatim
    phaseTurbulenceStabilisation
    {
        type    phaseTurbulenceStabilisation;

        libs    ("libmultiphaseEulerFoamFvModels.so");

        phase   air;

        alphaInversion  0.1;
    }
    \endverbatim

SourceFiles
    phaseTurbulenceStabilisation.C

\*---------------------------------------------------------------------------*/

#ifndef phaseTurbulenceStabilisation_H
#define phaseTurbulenceStabilisation_H

#include "fvModel.H"
#include "phaseModel.H"
#include "phaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class phaseTurbulenceStabilisation Declaration
\*---------------------------------------------------------------------------*/

class phaseTurbulenceStabilisation
:
    public fvModel
{
    // Private Data

        //- The name of the Lagrangian phase
        word phaseName_;

        //- Field names
        wordList fieldNames_;

        //- Phase-fraction below which turbulence property blending is applied
        dimensionedScalar alphaInversion_;

        //- Reference to the phase
        const phaseModel& phase_;

        //- Reference to the mixture turbulence model
        const phaseCompressible::momentumTransportModel& turbulence_;


    // Private Member Functions

        //- Add contribution to phase psi equation
        void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            tmp<volScalarField>
            (phaseCompressible::momentumTransportModel::*psi)() const
        ) const;


public:

    //- Runtime type information
    TypeName("phaseTurbulenceStabilisation");


    // Constructors

        //- Construct from explicit source name and mesh
        phaseTurbulenceStabilisation
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        phaseTurbulenceStabilisation
        (
            const phaseTurbulenceStabilisation&
        ) = delete;


    // Member Functions

        //- Return the list of fields for which the option adds source term
        //  to the transport equation
        virtual wordList addSupFields() const;

        //- Add contribution to phase k, epsilon or omega equation
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
        void operator=(const phaseTurbulenceStabilisation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
