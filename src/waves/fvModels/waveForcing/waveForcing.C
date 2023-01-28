/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "waveForcing.H"
#include "levelSet.H"
#include "fvMatrix.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

#include "fvcDdt.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(waveForcing, 0);
    addToRunTimeSelectionTable(fvModel, waveForcing, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::waveForcing::waveForcing
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    forcing(name, modelType, mesh, dict),
    waves_(waveSuperposition::New(mesh)),
    liquidPhaseName_(coeffs().lookup<word>("liquidPhase")),
    alphaName_(IOobject::groupName("alpha", liquidPhaseName_)),
    UName_(coeffs().lookupOrDefault<word>("U", "U")),
    scale_(this->scale())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::waveForcing::addSupFields() const
{
    return {alphaName_, UName_};
}


void Foam::fv::waveForcing::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (fieldName == alphaName_)
    {
        const volScalarField::Internal forceCoeff(this->forceCoeff(scale_));

        eqn -= fvm::Sp(forceCoeff, eqn.psi());
        eqn += forceCoeff*alphaWaves_();
    }
}


void Foam::fv::waveForcing::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (fieldName == UName_)
    {
        const volScalarField::Internal forceCoeff(rho*this->forceCoeff(scale_));

        eqn -= fvm::Sp(forceCoeff, eqn.psi());
        eqn += forceCoeff*Uwaves_();

        const surfaceScalarField& rhoPhi =
            mesh().lookupObject<surfaceScalarField>("rhoPhi");

        eqn += fvm::Sp
        (
            scale()*(fvc::ddt(rho)()() + fvc::div(rhoPhi)()()),
            eqn.psi()
        );
    }
}


bool Foam::fv::waveForcing::movePoints()
{
    scale_ = this->scale();
    return true;
}


void Foam::fv::waveForcing::topoChange(const polyTopoChangeMap&)
{
    scale_ = this->scale();
}


void Foam::fv::waveForcing::mapMesh(const polyMeshMap& map)
{
    scale_ = this->scale();
}


void Foam::fv::waveForcing::distribute(const polyDistributionMap&)
{
    scale_ = this->scale();
}


void Foam::fv::waveForcing::correct()
{
    const scalar t = mesh().time().value();

    // Cell centres and points
    const pointField& ccs = mesh().cellCentres();
    const pointField& pts = mesh().points();

    const scalarField h(waves_.height(t, ccs));
    const scalarField hp(waves_.height(t, pts));
    const vectorField uGas(waves_.UGas(t, ccs));
    const vectorField uGasp(waves_.UGas(t, pts));
    const vectorField uLiq(waves_.ULiquid(t, ccs));
    const vectorField uLiqp(waves_.ULiquid(t, pts));

    alphaWaves_ = volScalarField::Internal::New
    (
        "alphaWaves",
        mesh(),
        dimless,
        levelSetFraction(mesh(), h, hp, false)
    );

    Uwaves_ = volVectorField::Internal::New
    (
        "Uwaves",
        mesh(),
        dimVelocity,
        levelSetAverage(mesh(), h, hp, uGas, uGasp, uLiq, uLiqp)
    );
}


// ************************************************************************* //
