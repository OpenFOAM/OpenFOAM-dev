/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
#include "localEulerDdtScheme.H"
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
        autoPtr<fluidThermo>(psiuMulticomponentThermo::New(mesh).ptr())
    ),

    thermo_(refCast<psiuMulticomponentThermo>(isothermalFluid::thermo_)),

    b_(thermo_.Y("b")),

    unstrainedLaminarFlameSpeed(laminarFlameSpeed::New(thermo_)),

    Su
    (
        IOobject
        (
            "Su",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    SuMin(0.01*Su.average()),
    SuMax(4*Su.average()),

    Xi_
    (
        IOobject
        (
            "Xi",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    St
    (
        IOobject
        (
            "St",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Xi_*Su
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

    SuModel
    (
        combustionProperties.lookup("SuModel")
    ),

    sigmaExt("sigmaExt", dimless/dimTime, combustionProperties),

    XiModel
    (
        combustionProperties.lookup("XiModel")
    ),

    XiCoef("XiCoef", dimless, combustionProperties),

    XiShapeCoef("XiShapeCoef", dimless, combustionProperties),

    uPrimeCoef("uPrimeCoef", dimless, combustionProperties),

    ign(combustionProperties, runTime, mesh),

    thermophysicalTransport
    (
        momentumTransport(),
        thermo_,
        true
    ),

    thermo(thermo_),
    b(b_),
    Xi(Xi_)
{
    thermo.validate(type(), "ha", "ea");

    if (thermo_.containsSpecie("ft"))
    {
        fields.add(thermo_.Y("ft"));
    }

    fields.add(b);
    fields.add(thermo.he());
    fields.add(thermo.heu());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::XiFluid::~XiFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::XiFluid::prePredictor()
{
    isothermalFluid::prePredictor();

    if (pimple.predictTransport())
    {
        thermophysicalTransport.predict();
    }
}


void Foam::solvers::XiFluid::postCorrector()
{
    isothermalFluid::postCorrector();

    if (pimple.correctTransport())
    {
        thermophysicalTransport.correct();
    }
}


// ************************************************************************* //
