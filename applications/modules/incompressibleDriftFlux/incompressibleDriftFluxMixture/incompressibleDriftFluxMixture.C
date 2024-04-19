/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2024 OpenFOAM Foundation
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

#include "incompressibleDriftFluxMixture.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleDriftFluxMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleDriftFluxMixture::incompressibleDriftFluxMixture
(
    const fvMesh& mesh
)
:
    twoPhaseVoFMixture(mesh),

    nucModel_(viscosityModel::New(mesh, phase2Name())),
    muModel_(mixtureViscosityModel::New(*this)),

    rhod_("rho", dimDensity, muModel_()),
    rhoc_("rho", dimDensity, nucModel_()),

    alphaMax_(lookupOrDefault("alphaMax", 1.0)),

    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
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
    ),

    Uptr_(nullptr)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::incompressibleDriftFluxMixture::correct()
{
    rho_ = alpha1()*rhod_ + alpha2()*rhoc_;
    nu_ = muModel_->mu(rhoc_*nucModel_->nu(), *Uptr_)/rho_;
}


Foam::incompressibleDriftFluxMixture&
Foam::incompressibleDriftFluxMixture::initialise(const volVectorField& U)
{
    Uptr_ = &U;
    correct();
    return *this;
}


bool Foam::incompressibleDriftFluxMixture::read()
{
    if (twoPhaseMixture::read())
    {
        if (muModel_->read() || nucModel_->read())
        {
            nucModel_->lookup("rho") >> rhoc_;
            muModel_->lookup("rho") >> rhod_;

            alphaMax_ = muModel_->lookupOrDefault("alphaMax", 1.0);

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
