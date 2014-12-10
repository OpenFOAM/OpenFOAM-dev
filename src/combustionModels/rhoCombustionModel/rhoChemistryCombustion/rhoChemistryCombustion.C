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

#include "rhoChemistryCombustion.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::rhoChemistryCombustion::rhoChemistryCombustion
(
    const word& modelType,
    const fvMesh& mesh
)
:
    rhoCombustionModel(modelType, mesh),
    chemistryPtr_(rhoChemistryModel::New(mesh))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::rhoChemistryCombustion::~rhoChemistryCombustion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::rhoReactionThermo&
Foam::combustionModels::rhoChemistryCombustion::thermo()
{
    return chemistryPtr_->thermo();
}


const Foam::rhoReactionThermo&
Foam::combustionModels::rhoChemistryCombustion::thermo() const
{
    return chemistryPtr_->thermo();
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::rhoChemistryCombustion::rho() const
{
    return chemistryPtr_->thermo().rho();
}


// ************************************************************************* //
