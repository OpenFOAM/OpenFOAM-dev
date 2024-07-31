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

#include "fvTotalSource.H"
#include "fvCellSet.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvTotalSource, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvTotalSource::readCoeffs(const dictionary& dict)
{
    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fvTotalSource::addSource(fvMatrix<scalar>& eqn) const
{
    DebugInFunction
        << "eqnField=" << eqn.psi().name() << endl;

    const labelUList cells = this->cells();
    const scalar V = this->V();
    const dimensionedScalar S = this->S();

    // Check the dimensions
    eqn.dimensions() = S.dimensions();

    // Apply the source
    scalarField& eqnSource = eqn.source();
    forAll(cells, i)
    {
        const scalar f = mesh().V()[cells[i]]/V;
        eqnSource[cells[i]] -= f*S.value();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvTotalSource::fvTotalSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvSource(name, modelType, mesh, dict),
    phaseName_()
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvTotalSource::~fvTotalSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvTotalSource::addsSupToField(const word& fieldName) const
{
    const word group = IOobject::group(fieldName);

    return group == word::null || group == phaseName_;
}


Foam::tmp<Foam::scalarField> Foam::fvTotalSource::source
(
    const word& fieldName
) const
{
    return
        tmp<scalarField>
        (
            new scalarField
            (
                nCells(),
                (IOobject::group(fieldName) == phaseName_ ? S().value()/V() : 0)
            )
        );
}


bool Foam::fvTotalSource::read(const dictionary& dict)
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
