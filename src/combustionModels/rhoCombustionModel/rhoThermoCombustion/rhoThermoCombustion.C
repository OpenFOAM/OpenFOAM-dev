/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "rhoThermoCombustion.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::rhoThermoCombustion::rhoThermoCombustion
(
    const word& modelType,
    const fvMesh& mesh
)
:
    rhoCombustionModel(modelType, mesh),
    thermoPtr_(rhoReactionThermo::New(mesh))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::rhoThermoCombustion::~rhoThermoCombustion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::rhoReactionThermo&
Foam::combustionModels::rhoThermoCombustion::thermo()
{
    return thermoPtr_();
}


const Foam::rhoReactionThermo&
Foam::combustionModels::rhoThermoCombustion::thermo() const
{
    return thermoPtr_();
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::rhoThermoCombustion::rho() const
{
    return thermoPtr_().rho();
}


// ************************************************************************* //
