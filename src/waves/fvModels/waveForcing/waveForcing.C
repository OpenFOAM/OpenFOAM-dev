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


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::fv::waveForcing::readCoeffs()
{
    if (coeffs().found("lambdaCoeff"))
    {
        lambdaCoeff_ = coeffs().lookup<scalar>("lambdaCoeff");

        lambdaBoundaryCoeff_ =
            coeffs().lookupOrDefault<scalar>("lambdaBoundaryCoeff", 0);

        regionLength_ = this->regionLength();

        Info<< "    Volume average region length "
            << regionLength_.value() << endl;
    }
    else
    {
        readLambda();
    }
}


Foam::tmp<Foam::volScalarField::Internal> Foam::fv::waveForcing::scale() const
{
    return tmp<volScalarField::Internal>(scale_());
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::waveForcing::forceCoeff() const
{
    if (lambdaCoeff_ > 0)
    {
        const dimensionedScalar waveSpeed
        (
            dimVelocity,
            waves_.maxWaveSpeed(mesh().time().deltaTValue())
        );

        lambda_ = lambdaCoeff_*waveSpeed/regionLength_;
        lambdaBoundary_ = lambdaBoundaryCoeff_*waveSpeed/regionLength_;
    }

    return forcing::forceCoeff();
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
    lambdaCoeff_(0),
    lambdaBoundaryCoeff_(0),
    regionLength_("regionLength", dimLength, 0),
    waves_(waveSuperposition::New(mesh)),
    liquidPhaseName_(coeffs().lookup<word>("liquidPhase")),
    alphaName_(IOobject::groupName("alpha", liquidPhaseName_)),
    UName_(coeffs().lookupOrDefault<word>("U", "U")),
    scale_(forcing::scale().ptr())
{
    readCoeffs();
    writeForceFields();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::waveForcing::addSupFields() const
{
    return {alphaName_, UName_};
}


void Foam::fv::waveForcing::addSup
(
    const volScalarField& alpha,
    fvMatrix<scalar>& eqn
) const
{
    if (alpha.name() == alphaName_ && &alpha == &eqn.psi())
    {
        const volScalarField::Internal forceCoeff(this->forceCoeff());

        eqn -= fvm::Sp(forceCoeff, eqn.psi());
        eqn += forceCoeff*alphaWaves_();
    }
}


void Foam::fv::waveForcing::addSup
(
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    if (U.name() == UName_ && &U == &eqn.psi())
    {
        const volScalarField::Internal forceCoeff(rho*this->forceCoeff());

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
    // It is unlikely that the wave forcing region is moving
    // so this update could be removed or made optional
    scale_ = forcing::scale().ptr();
    return true;
}


void Foam::fv::waveForcing::topoChange(const polyTopoChangeMap&)
{
    scale_ = forcing::scale().ptr();
}


void Foam::fv::waveForcing::mapMesh(const polyMeshMap& map)
{
    scale_ = forcing::scale().ptr();
}


void Foam::fv::waveForcing::distribute(const polyDistributionMap&)
{
    scale_ = forcing::scale().ptr();
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::waveForcing::read(const dictionary& dict)
{
    if (forcing::read(dict))
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
