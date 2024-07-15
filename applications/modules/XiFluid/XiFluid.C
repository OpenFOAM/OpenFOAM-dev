/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

    SuModel_
    (
        SuModel::New
        (
            combustionProperties,
            thermo_,
            thermophysicalTransport
        )
    ),

    XiModel_
    (
        XiModel::New
        (
            combustionProperties,
            thermo_,
            thermophysicalTransport,
            SuModel_->Su()
        )
    ),

    ign(combustionProperties, runTime, mesh),

    thermo(thermo_),
    b(b_),
    Su(SuModel_->Su()),
    Xi(XiModel_->Xi())
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
