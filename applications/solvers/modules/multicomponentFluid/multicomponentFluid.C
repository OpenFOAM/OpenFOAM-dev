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

\*---------------------------------------------------------------------------*/

#include "multicomponentFluid.H"
#include "localEulerDdtScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(multicomponentFluid, 0);
    addToRunTimeSelectionTable(solver, multicomponentFluid, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::multicomponentFluid::multicomponentFluid(fvMesh& mesh)
:
    isothermalFluid
    (
        mesh,
        autoPtr<fluidThermo>(fluidMulticomponentThermo::New(mesh).ptr())
    ),

    thermo(refCast<fluidMulticomponentThermo>(isothermalFluid::thermo)),

    composition(thermo.composition()),

    Y(composition.Y()),

    reaction(combustionModel::New(thermo, turbulence())),

    thermophysicalTransport
    (
        fluidMulticomponentThermophysicalTransportModel::New
        (
            turbulence(),
            thermo
        )
    )
{
    forAll(Y, i)
    {
        fields.add(Y[i]);
    }
    fields.add(thermo.he());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::multicomponentFluid::~multicomponentFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::multicomponentFluid::thermophysicalTransportCorrector()
{
    if (pimple.transportCorr())
    {
        thermophysicalTransport->correct();
    }
}


// ************************************************************************* //
