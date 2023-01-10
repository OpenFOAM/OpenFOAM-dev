/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
        dimensionedScalar(dimViscosity, 0),
        calculatedFvPatchScalarField::typeName
    )
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::incompressibleMultiphaseVoFMixture::mu() const
{
    tmp<volScalarField> tmu
    (
        phases_[0]*phases_[0].rho()*phases_[0].nu()
    );
    volScalarField& mu = tmu.ref();

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        mu += phases_[phasei]*phases_[phasei].rho()*phases_[phasei].nu();
    }

    return tmu;
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleMultiphaseVoFMixture::mu(const label patchi) const
{
    tmp<scalarField> tmu
    (
        phases_[0].boundaryField()[patchi]
       *phases_[0].rho().value()
       *phases_[0].nu(patchi)
    );
    scalarField& mu = tmu.ref();

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        mu +=
            phases_[phasei].boundaryField()[patchi]
           *phases_[phasei].rho().value()
           *phases_[phasei].nu(patchi);
    }

    return tmu;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleMultiphaseVoFMixture::muf() const
{
    tmp<surfaceScalarField> tmuf
    (
        fvc::interpolate(phases_[0])
       *phases_[0].rho()*fvc::interpolate(phases_[0].nu())
    );
    surfaceScalarField& muf = tmuf.ref();

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        muf +=
            fvc::interpolate(phases_[phasei])
           *phases_[phasei].rho()*fvc::interpolate(phases_[phasei].nu());
    }

    return tmuf;
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleMultiphaseVoFMixture::nu() const
{
    return nu_;
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleMultiphaseVoFMixture::nu(const label patchi) const
{
    return nu_.boundaryField()[patchi];
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleMultiphaseVoFMixture::nuf() const
{
    return muf()/fvc::interpolate(rho());
}


void Foam::incompressibleMultiphaseVoFMixture::correct()
{
    forAll(phases_, phasei)
    {
        phases_[phasei].correct();
    }

    rho_ = phases_[0]*phases_[0].rho();

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        rho_ += phases_[phasei]*phases_[phasei].rho();
    }

    // Update the mixture kinematic viscosity
    nu_ = mu()/rho();

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
