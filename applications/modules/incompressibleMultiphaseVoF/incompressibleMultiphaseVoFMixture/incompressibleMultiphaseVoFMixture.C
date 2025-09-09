/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "incompressibleMultiphaseVoFMixture.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleMultiphaseVoFMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleMultiphaseVoFMixture::incompressibleMultiphaseVoFMixture
(
    const fvMesh& mesh
)
:
    multiphaseVoFMixture(mesh, incompressibleVoFphase::iNew(mesh)),

    phases_(multiphaseVoFMixture::phases().convert<incompressibleVoFphase>()),

    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 0)
    ),

    nu_
    (
        IOobject
        (
            "nu",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimKinematicViscosity, 0),
        calculatedFvPatchScalarField::typeName
    )
{
    correct();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::incompressibleMultiphaseVoFMixture::correct()
{
    forAll(phases_, phasei)
    {
        phases_[phasei].correct();
    }

    rho_ = phases_[0]*phases_[0].rho();
    volScalarField mu(rho_*phases_[0].nu());

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        const volScalarField alphaRho(phases_[phasei]*phases_[phasei].rho());
        rho_ += alphaRho;
        mu += alphaRho*phases_[phasei].nu();
    }

    // Update the mixture kinematic viscosity
    nu_ = mu/rho_;

    calcAlphas();
}


bool Foam::incompressibleMultiphaseVoFMixture::read()
{
    if (regIOobject::read())
    {
        lookup("sigmas") >> sigmas_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
