/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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

#include "KochFriedlanderSintering.H"
#include "fvmSup.H"
#include "fractal.H"
#include "addToRunTimeSelectionTable.H"
#include "sizeGroup.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(KochFriedlanderSintering, 0);
    addToRunTimeSelectionTable(fvModel, KochFriedlanderSintering, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dimensionSet Foam::fv::KochFriedlanderSintering::CsDims() const
{
    return pow(dimLength, -n_)*pow(dimTemperature, -m_)*dimTime;
}


void Foam::fv::KochFriedlanderSintering::readCoeffs(const dictionary& dict)
{
    if (dict.lookup<word>("populationBalance") != popBal_.name())
    {
        FatalIOErrorInFunction(dict)
            << "Cannot change the population balance of a " << type()
            << " model at run time" << exit(FatalIOError);
    }

    n_ = dict.lookup<scalar>("n");
    m_ = dict.lookup<scalar>("m");
    Cs_.read(dict, CsDims());
    Ta_.read(dict);
    dpMin_.readIfPresent(dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::KochFriedlanderSintering::KochFriedlanderSintering
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    popBal_
    (
        mesh().lookupObject<diameterModels::populationBalanceModel>
        (
            coeffs(dict).lookup("populationBalance")
        )
    ),
    n_(coeffs(dict).lookup<scalar>("n")),
    m_(coeffs(dict).lookup<scalar>("m")),
    Cs_("Cs", CsDims(), coeffs(dict)),
    Ta_("Ta", dimTemperature, coeffs(dict)),
    dpMin_("dpMin", dimLength, coeffs(dict), scalar(0))
{
    readCoeffs(coeffs(dict));

    forAll(popBal_.sizeGroups(), sizeGroupi)
    {
        const diameterModels::shapeModel& shapeModel =
            popBal_.sizeGroups()[sizeGroupi].shape();

        if (!isA<diameterModels::shapeModels::fractal>(shapeModel)) continue;

        const diameterModels::shapeModels::fractal& fractal =
            refCast<const diameterModels::shapeModels::fractal>(shapeModel);

        kappaNameToSizeGroupIndices_.insert(fractal.fld().name(), sizeGroupi);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::KochFriedlanderSintering::addsSupToField
(
    const word& fieldName
) const
{
    return kappaNameToSizeGroupIndices_.found(fieldName);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::KochFriedlanderSintering::tau
(
    const volScalarField::Internal& kappa
) const
{
    const diameterModels::sizeGroup& fi =
        popBal_.sizeGroups()[kappaNameToSizeGroupIndices_[kappa.name()]];

    const volScalarField::Internal& T = fi.phase().thermo().T();

    const volScalarField::Internal dp(6/max(6/fi.dSph(), kappa));

    return Cs_*pow(dp, n_)*pow(T, m_)*exp(Ta_/T*(1 - dpMin_/dp));
}


void Foam::fv::KochFriedlanderSintering::addSup
(
    const volScalarField& alphaFi,
    const volScalarField& rho,
    const volScalarField& kappa,
    fvMatrix<scalar>& eqn
) const
{
    const diameterModels::sizeGroup& fi =
        popBal_.sizeGroups()[kappaNameToSizeGroupIndices_[kappa.name()]];

    const volScalarField::Internal R(alphaFi()*rho()/tau(kappa));

    eqn += R*6/fi.dSph() - fvm::Sp(R, kappa);
}


bool Foam::fv::KochFriedlanderSintering::movePoints()
{
    return true;
}


void Foam::fv::KochFriedlanderSintering::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::KochFriedlanderSintering::mapMesh(const polyMeshMap&)
{}


void Foam::fv::KochFriedlanderSintering::distribute(const polyDistributionMap&)
{}


bool Foam::fv::KochFriedlanderSintering::read(const dictionary& dict)
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
