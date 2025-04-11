/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2025 OpenFOAM Foundation
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

void Foam::fv::coefficientPhaseChange::readCoeffs(const dictionary& dict)
{
    reReadSpecies(dict);

    C_.read(dict);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::coefficientPhaseChange::timesY1
(
    tmp<volScalarField::Internal> mDot
) const
{
    const ThermoRefPair<multicomponentThermo> mcThermos =
        thermos().thermos<multicomponentThermo>();

    if (!mcThermos.valid().first() || species().empty())
    {
        return mDot;
    }

    if (species().size() == 1)
    {
        return mcThermos.first().Y()[specieis().first()]*mDot;
    }

    tmp<volScalarField::Internal> tY1 =
        volScalarField::Internal::New
        (
            typedName("Y1"),
            mcThermos.first().Y()[specieis(0).first()]
        );

    for (label mDoti = 1; mDoti < species().size(); ++ mDoti)
    {
        tY1.ref() += mcThermos.first().Y()[specieis(mDoti).first()];
    }

    return tY1*mDot;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::coefficientPhaseChange::mDotByAlpha1Y1() const
{
    return C_*mag(fvc::grad(alpha1_))()();
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
    phaseChange
    (
        name,
        modelType,
        mesh,
        dict,
        readSpecies(coeffs(modelType, dict), false)
    ),
    C_("C", dimMass/dimArea/dimTime, NaN),
    alpha1_(mesh().lookupObject<volScalarField>(alphaNames().first()))
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::coefficientPhaseChange::mDot() const
{
    return timesY1(alpha1_*mDotByAlpha1Y1());
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::coefficientPhaseChange::mDot(const label mDoti) const
{
    const ThermoRefPair<multicomponentThermo> mcThermos =
        thermos().thermos<multicomponentThermo>();

    tmp<volScalarField::Internal> tmDot = alpha1_*mDotByAlpha1Y1();

    if (mcThermos.valid().first())
    {
        const labelPair specieis = this->specieis(mDoti);

        tmDot.ref() *= mcThermos.first().Y()[specieis.first()];
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
        const volScalarField::Internal mDotByAlpha1(timesY1(mDotByAlpha1Y1()));

        if (i == 0)
        {
            eqn -= fvm::Sp(mDotByAlpha1, eqn.psi());
        }
        else
        {
            eqn +=
                mDotByAlpha1*alpha1_
              - correction(fvm::Sp(mDotByAlpha1, eqn.psi()));
        }
    }
    else
    {
        massTransfer::addSup(alpha, rho, eqn);
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
    const label i = index(phaseNames(), eqn.psi().group());

    const ThermoRefPair<multicomponentThermo>& mcThermos =
        thermos().thermos<multicomponentThermo>();

    if
    (
        i == 0
     && mcThermos.valid().first()
     && mcThermos.first().containsSpecie(Yi.member())
     && species().found(Yi.member())
    )
    {
        eqn -= fvm::Sp(alpha1_*mDotByAlpha1Y1(), Yi);
    }
    else
    {
        phaseChange::addSup(alpha, rho, Yi, eqn);
    }
}


bool Foam::fv::coefficientPhaseChange::read(const dictionary& dict)
{
    if (phaseChange::read(dict))
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
