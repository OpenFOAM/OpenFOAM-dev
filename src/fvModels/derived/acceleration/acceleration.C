/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2024 OpenFOAM Foundation
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
#include "acceleration.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(acceleration, 0);
    addToRunTimeSelectionTable(fvModel, acceleration, dictionary);
    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvModel,
        acceleration,
        dictionary,
        accelerationSource,
        "accelerationSource"
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::acceleration::readCoeffs()
{
    UName_ = coeffs().lookupOrDefault<word>("U", "U");

    velocity_ = Function1<vector>::New("velocity", coeffs());
}


template<class AlphaRhoFieldType>
void Foam::fv::acceleration::add
(
    const AlphaRhoFieldType& alphaRho,
    fvMatrix<vector>& eqn
) const
{
    const DimensionedField<scalar, volMesh>& V = mesh().V();

    const scalar t = mesh().time().value();
    const scalar dt = mesh().time().deltaTValue();
    const vector dU =
        velocity_->value(mesh().time().timeToUserTime(t))
      - velocity_->value(mesh().time().timeToUserTime(t - dt));
    const vector a = dU/mesh().time().deltaTValue();

    const labelUList cells = set_.cells();

    vectorField& eqnSource = eqn.source();
    forAll(cells, i)
    {
        eqnSource[cells[i]] -= V[cells[i]]*alphaRho[cells[i]]*a;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::acceleration::acceleration
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    set_(mesh, coeffs()),
    UName_(word::null),
    velocity_(nullptr)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::acceleration::addSupFields() const
{
    return wordList(1, UName_);
}


void Foam::fv::acceleration::addSup
(
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    add(geometricOneField(), eqn);
}


void Foam::fv::acceleration::addSup
(
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    add(rho, eqn);
}


void Foam::fv::acceleration::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    add((alpha*rho)(), eqn);
}


bool Foam::fv::acceleration::movePoints()
{
    set_.movePoints();
    return true;
}


void Foam::fv::acceleration::topoChange(const polyTopoChangeMap& map)
{
    set_.topoChange(map);
}


void Foam::fv::acceleration::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void Foam::fv::acceleration::distribute(const polyDistributionMap& map)
{
    set_.distribute(map);
}


bool Foam::fv::acceleration::read(const dictionary& dict)
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
