/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "PaSR.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(PaSR, 0);
    addToRunTimeSelectionTable(combustionModel, PaSR, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::PaSR::PaSR
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
:
    laminar(modelType, thermo, turb, combustionProperties),
    Cmix_(this->coeffs().template lookup<scalar>("Cmix")),
    kappa_
    (
        IOobject
        (
            thermo.phasePropertyName(typeName + ":kappa"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::PaSR::~PaSR()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::PaSR::correct()
{
    laminar::correct();

    tmp<volScalarField> tepsilon(this->turbulence().epsilon());
    const scalarField& epsilon = tepsilon();

    tmp<volScalarField> tnuEff(this->turbulence().nuEff());
    const scalarField& nuEff = tnuEff();

    tmp<volScalarField> ttc(this->chemistryPtr_->tc());
    const scalarField& tc = ttc();

    forAll(epsilon, i)
    {
        const scalar tk =
            Cmix_*sqrt(max(nuEff[i]/(epsilon[i] + small), 0));

        if (tk > small)
        {
            kappa_[i] = tc[i]/(tc[i] + tk);
        }
        else
        {
            kappa_[i] = 1.0;
        }
    }
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::PaSR::R(volScalarField& Y) const
{
    return kappa_*laminar::R(Y);
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::PaSR::Qdot() const
{
    return volScalarField::New
    (
        this->thermo().phasePropertyName(typeName + ":Qdot"),
        kappa_*laminar::Qdot()
    );
}


bool Foam::combustionModels::PaSR::read()
{
    if (laminar::read())
    {
        this->coeffs().lookup("Cmix") >> Cmix_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
