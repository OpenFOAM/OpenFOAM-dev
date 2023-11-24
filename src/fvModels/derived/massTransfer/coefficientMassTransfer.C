/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "coefficientMassTransfer.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(coefficientMassTransfer, 0);
    addToRunTimeSelectionTable(fvModel, coefficientMassTransfer, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::coefficientMassTransfer::readCoeffs()
{
    C_.read(coeffs());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::coefficientMassTransfer::coefficientMassTransfer
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    massTransferBase(name, modelType, mesh, dict),
    C_("C", dimMass/dimArea/dimTime, NaN)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::coefficientMassTransfer::mDot() const
{
    const volScalarField& alpha1 =
        mesh().lookupObject<volScalarField>(alphaNames().first());

    return C_*alpha1()*mag(fvc::grad(alpha1))()();
}


void Foam::fv::coefficientMassTransfer::addSup
(
    const volScalarField& alpha,
    fvMatrix<scalar>& eqn
) const
{
    const label i = index(alphaNames(), eqn.psi().name());

    if (i != -1)
    {
        const volScalarField& alpha1 =
            mesh().lookupObject<volScalarField>(alphaNames().first());

        const volScalarField::Internal SByAlpha1
        (
            C_*mag(fvc::grad(alpha1)()())/rho(i)
        );

        if (i == 0)
        {
            eqn -= fvm::Sp(SByAlpha1, eqn.psi());
        }
        else
        {
            eqn += SByAlpha1*alpha1 - correction(fvm::Sp(SByAlpha1, eqn.psi()));
        }
    }
    else
    {
        massTransferBase::addSupType(alpha, eqn);
    }
}


void Foam::fv::coefficientMassTransfer::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    const label i = index(alphaNames(), eqn.psi().name());

    if (i != -1)
    {
        const volScalarField& alpha1 =
            mesh().lookupObject<volScalarField>(alphaNames().first());

        const volScalarField::Internal SByAlpha1
        (
            C_*mag(fvc::grad(alpha1))
        );

        if (i == 0)
        {
            eqn -= fvm::Sp(SByAlpha1, eqn.psi());
        }
        else
        {
            eqn += SByAlpha1*alpha1 - correction(fvm::Sp(SByAlpha1, eqn.psi()));
        }
    }
    else
    {
        massTransferBase::addSupType(alpha, rho, eqn);
    }
}


bool Foam::fv::coefficientMassTransfer::read(const dictionary& dict)
{
    if (massTransferBase::read(dict))
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
