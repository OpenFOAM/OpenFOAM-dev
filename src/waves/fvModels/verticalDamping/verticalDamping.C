/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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

#include "verticalDamping.H"
#include "fvMatrix.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(verticalDamping, 0);
    addToRunTimeSelectionTable(fvModel, verticalDamping, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::verticalDamping::readCoeffs()
{
    readLambda();
}


void Foam::fv::verticalDamping::add
(
    const volVectorField& alphaRhoU,
    fvMatrix<vector>& eqn
) const
{
    const uniformDimensionedVectorField& g =
        mesh().lookupObject<uniformDimensionedVectorField>("g");

    const dimensionedSymmTensor gg(sqr(g)/magSqr(g));

    eqn -= forceCoeff()*(gg & alphaRhoU());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::verticalDamping::verticalDamping
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    forcing(name, modelType, mesh, dict),
    UName_(coeffs().lookupOrDefault<word>("U", "U"))
{
    writeForceFields();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::verticalDamping::addSupFields() const
{
    return wordList(1, UName_);
}


void Foam::fv::verticalDamping::addSup
(
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    add(U, eqn);
}


void Foam::fv::verticalDamping::addSup
(
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    add(rho*U, eqn);
}


void Foam::fv::verticalDamping::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    add(alpha*rho*U, eqn);
}


bool Foam::fv::verticalDamping::movePoints()
{
    return true;
}


void Foam::fv::verticalDamping::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::verticalDamping::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::verticalDamping::distribute(const polyDistributionMap&)
{}


bool Foam::fv::verticalDamping::read(const dictionary& dict)
{
    if (forcing::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
