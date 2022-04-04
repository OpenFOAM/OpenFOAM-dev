/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2022 OpenFOAM Foundation
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

#include "isotropicDamping.H"
#include "fvMatrix.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(isotropicDamping, 0);
    addToRunTimeSelectionTable(fvModel, isotropicDamping, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::isotropicDamping::readCoeffs()
{
    value_ =
        dimensionedVector
        (
            value_.name(),
            value_.dimensions(),
            coeffs().lookup(value_.name())
        );
}


void Foam::fv::isotropicDamping::add
(
    const volScalarField::Internal& forceCoeff,
    fvMatrix<vector>& eqn
) const
{
    eqn -= fvm::Sp(forceCoeff, eqn.psi());
    eqn += forceCoeff*value_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::isotropicDamping::isotropicDamping
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    damping(name, modelType, dict, mesh),
    value_("value", dimVelocity, vector::uniform(NaN))
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::isotropicDamping::addSup
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    add(this->forceCoeff(), eqn);
}


void Foam::fv::isotropicDamping::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    add(rho*forceCoeff(), eqn);
}


void Foam::fv::isotropicDamping::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    add(alpha()*rho()*this->forceCoeff(), eqn);
}


bool Foam::fv::isotropicDamping::movePoints()
{
    return true;
}


void Foam::fv::isotropicDamping::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::isotropicDamping::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::isotropicDamping::distribute(const polyDistributionMap&)
{}


bool Foam::fv::isotropicDamping::read(const dictionary& dict)
{
    if (damping::read(dict))
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
