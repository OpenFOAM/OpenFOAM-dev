/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
#include "fvMesh.H"
#include "fvMatrix.H"
#include "geometricOneField.H"
#include "meshTools.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(verticalDamping, 0);
    addToRunTimeSelectionTable(option, verticalDamping, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::verticalDamping::add
(
    const volVectorField& alphaRhoU,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const uniformDimensionedVectorField& g =
        mesh_.lookupObject<uniformDimensionedVectorField>("g");

    const dimensionedSymmTensor lgg(lambda_*sqr(g)/magSqr(g));

    const DimensionedField<scalar, volMesh>& V = mesh_.V();

    // Check dimensions
    eqn.dimensions()
      - V.dimensions()*(lgg.dimensions() & alphaRhoU.dimensions());

    // Calculate the force and apply it to the equation
    vectorField force(cells_.size());
    forAll(cells_, i)
    {
        const label c = cells_[i];
        force[i] = V[c]*(lgg.value() & alphaRhoU[c]);
    }
    meshTools::constrainDirection(mesh_, mesh_.solutionD(), force);
    forAll(cells_, i)
    {
        const label c = cells_[i];
        eqn.source()[c] += force[i];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::verticalDamping::verticalDamping
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    lambda_("lambda", dimless/dimTime, coeffs_.lookup("lambda"))
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::verticalDamping::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    add(eqn.psi(), eqn, fieldi);
}


void Foam::fv::verticalDamping::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    add(rho*eqn.psi(), eqn, fieldi);
}


void Foam::fv::verticalDamping::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    add(alpha*rho*eqn.psi(), eqn, fieldi);
}


bool Foam::fv::verticalDamping::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        lambda_ =
            dimensionedScalar
            (
                lambda_.name(),
                lambda_.dimensions(),
                coeffs_.lookup(lambda_.name())
            );

        fieldNames_ = wordList(1, coeffs_.lookupOrDefault<word>("U", "U"));

        applied_.setSize(1, false);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
