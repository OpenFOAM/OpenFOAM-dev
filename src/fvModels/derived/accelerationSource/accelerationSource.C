/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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
    addToRunTimeSelectionTable(fvModel, accelerationSource, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::accelerationSource::readCoeffs()
{
    UName_ = coeffs().lookupOrDefault<word>("U", "U");

    velocity_ = Function1<vector>::New("velocity", coeffs());
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
    fvModel(name, modelType, dict, mesh),
    set_(coeffs(), mesh),
    UName_(word::null),
    velocity_(nullptr)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::accelerationSource::addSupFields() const
{
    return wordList(1, UName_);
}


void Foam::fv::accelerationSource::addSup
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    add(geometricOneField(), eqn, fieldName);
}


void Foam::fv::accelerationSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    add(rho, eqn, fieldName);
}


void Foam::fv::accelerationSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    add((alpha*rho)(), eqn, fieldName);
}


void Foam::fv::accelerationSource::updateMesh(const mapPolyMesh& map)
{
    set_.updateMesh(map);
}


void Foam::fv::accelerationSource::distribute(const mapDistributePolyMesh& map)
{
    set_.distribute(map);
}


bool Foam::fv::accelerationSource::movePoints()
{
    set_.movePoints();
    return true;
}


bool Foam::fv::accelerationSource::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        set_.read(coeffs());
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
