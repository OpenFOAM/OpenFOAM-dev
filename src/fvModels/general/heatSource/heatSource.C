/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "heatSource.H"
#include "basicThermo.H"
#include "fvModels.H"
#include "fvMatrix.H"
#include "Scale.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(heatSource, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        heatSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::heatSource::readCoeffs(const dictionary& dict)
{
    if (!dict.found("q") && !dict.found("Q"))
    {
        FatalIOErrorInFunction(coeffs(dict))
            << "Neither heat source per unit volume, q, or total heat source, "
            << "Q, has been specified. One is required." << exit(FatalIOError);
    }

    if (dict.found("q") && dict.found("Q"))
    {
        FatalIOErrorInFunction(coeffs(dict))
            << "Both heat source per unit volume, q, and total heat source, "
            << "Q, have been specified. One is required."
            << exit(FatalIOError);
    }

    if (dict.found("q"))
    {
        q_.reset
        (
            Function1<scalar>::New
            (
                "q",
                mesh().time().userUnits(),
                dimPower/dimVolume,
                dict
            ).ptr()
        );
    }
    else
    {
        q_.reset
        (
            new Function1s::Scale<scalar>
            (
                "q",
                Function1s::Constant<scalar>("1/V", 1/set_.V()),
                Function1s::Constant<scalar>("1", 1),
                Function1<scalar>::New
                (
                    "Q",
                    mesh().time().userUnits(),
                    dimPower,
                    dict
                )()
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatSource::heatSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    set_(mesh, coeffs(dict)),
    q_(nullptr)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::heatSource::~heatSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::heatSource::addSupFields() const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>(physicalProperties::typeName);

    return wordList(1, thermo.he().name());
}


void Foam::fv::heatSource::addSup
(
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    const labelUList cells = set_.cells();

    const scalar q = q_->value(mesh().time().value());

    scalarField& eqnSource = eqn.source();
    forAll(cells, i)
    {
        eqnSource[cells[i]] -= mesh().V()[cells[i]]*q;
    }
}


void Foam::fv::heatSource::addSup
(
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    addSup(he, eqn);
}


bool Foam::fv::heatSource::movePoints()
{
    set_.movePoints();
    return true;
}


void Foam::fv::heatSource::topoChange(const polyTopoChangeMap& map)
{
    set_.topoChange(map);
}


void Foam::fv::heatSource::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void Foam::fv::heatSource::distribute(const polyDistributionMap& map)
{
    set_.distribute(map);
}


bool Foam::fv::heatSource::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
