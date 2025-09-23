/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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

#include "XiFluid.H"
#include "bXiIgnition.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(XiFluid, 0);
    addToRunTimeSelectionTable(solver, XiFluid, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::XiFluid::XiFluid(fvMesh& mesh)
:
    isothermalFluid
    (
        mesh,
        autoPtr<fluidThermo>(new ubPsiThermo(mesh))
    ),

    thermo_(refCast<ubPsiThermo>(isothermalFluid::thermo_)),

    thermophysicalTransport
    (
        momentumTransport(),
        thermo_,
        true
    ),

    combustionProperties
    (
        IOobject
        (
            "combustionProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    ignited_(false),

    bMin_
    (
        combustionProperties.subDict("flameWrinkling")
       .lookupOrDefault("bMin", 1e-3)
    ),

    mgbCoeff_
    (
        combustionProperties.subDict("flameWrinkling")
       .lookupOrDefault("mgbCoeff", 1e-3)
    ),

    SuModel_
    (
        SuModel::New
        (
            combustionProperties,
            thermo_.uThermo(),
            momentumTransport()
        )
    ),

    XiModel_
    (
        XiModel::New
        (
            combustionProperties,
            thermo_,
            momentumTransport(),
            SuModel_->Su()
        )
    ),

    thermo(thermo_),

    b(thermo.b()),
    uThermo(thermo.uThermo()),

    c(thermo.c()),
    bThermo(thermo.bThermo()),

    Su(SuModel_->Su()),
    Xi(XiModel_->Xi())
{
    mesh.schemes().setFluxRequired(b.name());

    const UPtrListDictionary<fv::bXiIgnition> ignitionModels
    (
        fvModels().lookupType<fv::bXiIgnition>()
    );

    forAll(ignitionModels, i)
    {
        if (ignitionModels[i].ignited())
        {
            ignited_ = true;
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::XiFluid::~XiFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::XiFluid::thermophysicalTransportPredictor()
{
    thermophysicalTransport.predict();
}


void Foam::solvers::XiFluid::thermophysicalTransportCorrector()
{
    thermophysicalTransport.correct();
}


void Foam::solvers::XiFluid::reset()
{
    ignited_ = false;
    thermo_.reset();

    const surfaceScalarField phib("phib", phi);
    thermo_.b().correctBoundaryConditions();
    thermo_.c() = 1.0 - thermo_.b();

    SuModel_->reset();
    XiModel_->reset();
}

// ************************************************************************* //
