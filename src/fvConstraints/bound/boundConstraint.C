/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "boundConstraint.H"
#include "volFields.H"
#include "bound.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(bound, 0);
    addToRunTimeSelectionTable
    (
        fvConstraint,
        bound,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::bound::readCoeffs()
{
    fieldName_ = coeffs().lookup<word>("field");
    min_ = coeffs().lookup<scalar>("min");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::bound::bound
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvConstraint(name, modelType, mesh, dict),
    fieldName_(word::null),
    min_(0)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::bound::constrainedFields() const
{
    return wordList(1, fieldName_);
}


bool Foam::fv::bound::constrain(volScalarField& f) const
{
    return Foam::bound(f, dimensionedScalar(f.dimensions(), min_));
}


bool Foam::fv::bound::movePoints()
{
    return true;
}


void Foam::fv::bound::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::bound::mapMesh(const polyMeshMap&)
{}


void Foam::fv::bound::distribute(const polyDistributionMap&)
{}


bool Foam::fv::bound::read(const dictionary& dict)
{
    if (fvConstraint::read(dict))
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
