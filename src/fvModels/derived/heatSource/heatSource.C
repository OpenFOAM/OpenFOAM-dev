/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

void Foam::fv::heatSource::readCoeffs()
{
    if (!coeffs().found("q") && !coeffs().found("Q"))
    {
        FatalIOErrorInFunction(coeffs())
            << "Neither heat source per unit volume, q, or total heat source, "
            << "Q, has been specified. One is required." << exit(FatalIOError);
    }

    if (coeffs().found("q") && coeffs().found("Q"))
    {
        FatalIOErrorInFunction(coeffs())
            << "Both heat source per unit volume, q, and total heat source, "
            << "Q, have been specified. One is required."
            << exit(FatalIOError);
    }

    if (coeffs().found("q"))
    {
        q_.reset(Function1<scalar>::New("q", coeffs()).ptr());
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
                Function1<scalar>::New("Q", coeffs())()
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatSource::heatSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    set_(coeffs(), mesh),
    q_(nullptr)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::heatSource::~heatSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::heatSource::addSupFields() const
{
    const basicThermo& thermo =
        mesh().lookupObject<basicThermo>(basicThermo::dictName);

    return wordList(1, thermo.he().name());
}


void Foam::fv::heatSource::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    const labelList& cells = set_.cells();

    const scalar t = mesh().time().value();
    const scalar q = q_->value(t);

    forAll(cells, i)
    {
        eqn.source()[cells[i]] -= mesh().V()[cells[i]]*q;
    }
}


void Foam::fv::heatSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    addSup(eqn, fieldName);
}


bool Foam::fv::heatSource::read(const dictionary& dict)
{
    if (fvModel::read(dict))
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
