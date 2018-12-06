/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "fvMesh.H"
#include "fvMatrix.H"
#include "geometricOneField.H"
#include "accelerationSource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(accelerationSource, 0);
    addToRunTimeSelectionTable(option, accelerationSource, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::accelerationSource::accelerationSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    velocity_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::accelerationSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    add(geometricOneField(), eqn, fieldi);
}


void Foam::fv::accelerationSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    add(rho, eqn, fieldi);
}


void Foam::fv::accelerationSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    add((alpha*rho)(), eqn, fieldi);
}


bool Foam::fv::accelerationSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        fieldNames_ = wordList(1, coeffs_.lookupOrDefault<word>("U", "U"));

        applied_.setSize(1, false);

        velocity_ = Function1<vector>::New("velocity", dict);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
