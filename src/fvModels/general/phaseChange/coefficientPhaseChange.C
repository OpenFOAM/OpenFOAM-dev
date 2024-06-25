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

#include "coefficientPhaseChange.H"
#include "fvcGrad.H"
#include "multicomponentThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(coefficientPhaseChange, 0);
    addToRunTimeSelectionTable(fvModel, coefficientPhaseChange, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::coefficientPhaseChange::readCoeffs()
{
    C_.read(coeffs());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::coefficientPhaseChange::coefficientPhaseChange
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    singleComponentPhaseChange
    (
        name,
        modelType,
        mesh,
        dict,
        {false, false},
        {false, false}
    ),
    C_("C", dimMass/dimArea/dimTime, NaN)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::coefficientPhaseChange::mDot() const
{
    const volScalarField& alpha1 =
        mesh().lookupObject<volScalarField>(alphaNames().first());

    tmp<volScalarField::Internal> tmDot =
        C_*alpha1()*mag(fvc::grad(alpha1))()();

    if (specieis().first() != -1)
    {
        tmDot.ref() *= specieThermos().first().Y()[specieis().first()];
    }

    return tmDot;
}


void Foam::fv::coefficientPhaseChange::addSup
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

        volScalarField::Internal mDotByAlpha1(C_*mag(fvc::grad(alpha1)));

        if (specieis().first() != -1)
        {
            mDotByAlpha1 *= specieThermos().first().Y()[specieis().first()];
        }

        if (i == 0)
        {
            eqn -= fvm::Sp(mDotByAlpha1, eqn.psi());
        }
        else
        {
            eqn +=
                mDotByAlpha1*alpha1
              - correction(fvm::Sp(mDotByAlpha1, eqn.psi()));
        }
    }
    else
    {
        phaseChange::addSup(alpha, rho, eqn);
    }
}


void Foam::fv::coefficientPhaseChange::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& Yi,
    fvMatrix<scalar>& eqn
) const
{
    const label i = index(alphaNames(), eqn.psi().name());

    if (i == 0 && specieis().first() != -1 && Yi.member() == specie())
    {
        eqn -= fvm::Sp(C_*alpha()*mag(fvc::grad(alpha))()(), Yi);
    }
    else
    {
        singleComponentPhaseChange::addSup(alpha, rho, Yi, eqn);
    }
}


bool Foam::fv::coefficientPhaseChange::read(const dictionary& dict)
{
    if (singleComponentPhaseChange::read(dict))
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
