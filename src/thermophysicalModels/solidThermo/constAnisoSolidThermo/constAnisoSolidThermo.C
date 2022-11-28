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

#include "constAnisoSolidThermo.H"
#include "addToRunTimeSelectionTable.H"

/* * * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(constAnisoSolidThermo, 0);
    addToRunTimeSelectionTable(basicThermo, constAnisoSolidThermo, fvMesh);
    addToRunTimeSelectionTable(solidThermo, constAnisoSolidThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constAnisoSolidThermo::constAnisoSolidThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    constSolidThermo(mesh, false, phaseName),
    Kappa_(readProperty<vector>("Kappa", kappa_.dimensions()))
{
    kappa_ = mag(Kappa_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constAnisoSolidThermo::~constAnisoSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::constAnisoSolidThermo::kappa() const
{
    NotImplemented;
    return kappa_;
}


const Foam::volVectorField& Foam::constAnisoSolidThermo::Kappa() const
{
    return Kappa_;
}


// ************************************************************************* //
