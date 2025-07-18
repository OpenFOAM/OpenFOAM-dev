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

#include "homogeneousCondensation.H"
#include "fundamentalConstants.H"
#include "multicomponentThermo.H"
#include "physicoChemicalConstants.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(homogeneousCondensation, 0);
    addToRunTimeSelectionTable(fvModel, homogeneousCondensation, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::homogeneousCondensation::readCoeffs(const dictionary& dict)
{
    reReadSpecie(dict);

    saturationModel_.reset
    (
        saturationPressureModel::New("pSat", dict).ptr()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::homogeneousCondensation::homogeneousCondensation
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phaseChange(name, modelType, mesh, dict, readSpecie(dict, true)),
    fluid_
    (
        mesh().lookupObject<phaseSystem>(phaseSystem::propertiesName)
    ),
    d_
    (
        IOobject
        (
            name + ":d",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimLength, 0)
    ),
    mDotByAlphaGas_
    (
        IOobject
        (
            name + ":mDotByAlpha",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimDensity/dimTime, 0)
    ),
    saturationModel_(nullptr)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::homogeneousCondensation::d() const
{
    return d_;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::homogeneousCondensation::nDot() const
{
    const volScalarField::Internal& alphaGas =
        mesh().lookupObject<volScalarField::Internal>(alphaNames().first());

    const ThermoRefPair<multicomponentThermo> multicomponentThermos =
        this->multicomponentThermos(true, false);

    const multicomponentThermo& thermoGas = multicomponentThermos.first();

    const volScalarField& p = this->p();
    const volScalarField& T = thermoGas.T();

    const volScalarField::Internal rhoLiquid
    (
        multicomponentThermos.valid().second()
      ? vfToVif(multicomponentThermos.second().rhoi(specieis().second(), p, T))
      : vfToVif(thermos().second().rho())
    );

    const volScalarField::Internal v(constant::mathematical::pi/6*pow3(d_));

    return alphaGas*mDotByAlphaGas_/(rhoLiquid*v);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::homogeneousCondensation::mDot() const
{
    const volScalarField::Internal& alphaGas =
        mesh().lookupObject<volScalarField::Internal>(alphaNames().first());

    return alphaGas*mDotByAlphaGas_;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::homogeneousCondensation::tau() const
{
    static const dimensionedScalar mDotRootVSmall
    (
        dimDensity/dimTime,
        rootVSmall
    );

    const ThermoRefPair<multicomponentThermo> multicomponentThermos =
        this->multicomponentThermos(true, false);

    const multicomponentThermo& thermoGas = multicomponentThermos.first();

    // Phase molecular masses and densities
    const volScalarField::Internal WGas(vfToVif(thermoGas.W()));
    const volScalarField::Internal rhoGas(vfToVif(thermoGas.rho()));

    // Mole fraction of nucleating specie
    const volScalarField::Internal Xi
    (
        thermoGas.Y()[specieis().first()]*WGas/thermoGas.Wi(specieis().first())
    );

    return Xi*rhoGas/max(mDotByAlphaGas_, mDotRootVSmall);
}


void Foam::fv::homogeneousCondensation::correct()
{
    #define infoFieldVariable(field, print) infoField(#field, field, print)

    Info<< type() << ": " << name() << endl << incrIndent;

    using constant::mathematical::pi;
    using constant::physicoChemical::NNA;
    using constant::physicoChemical::k;

    const ThermoRefPair<multicomponentThermo> multicomponentThermos =
        this->multicomponentThermos(true, false);

    const multicomponentThermo& thermoGas = multicomponentThermos.first();

    const volScalarField& p = this->p();
    const volScalarField& T = thermoGas.T();
    infoFieldVariable(p, debug);
    infoFieldVariable(T, debug);

    // Phase molecular masses and densities
    const volScalarField::Internal WGas(vfToVif(thermoGas.W()));
    const volScalarField::Internal rhoGas(vfToVif(thermoGas.rho()));
    const volScalarField::Internal WLiquid
    (
        multicomponentThermos.valid().second()
      ? volScalarField::Internal::New
        (
            "W",
            mesh(),
            multicomponentThermos.second().Wi(specieis().second())
        )
      : vfToVif(thermos().second().W())
    );
    const volScalarField::Internal rhoLiquid
    (
        multicomponentThermos.valid().second()
      ? vfToVif(multicomponentThermos.second().rhoi(specieis().second(), p, T))
      : vfToVif(thermos().second().rho())
    );
    infoFieldVariable(WGas, debug);
    infoFieldVariable(rhoGas, debug);
    infoFieldVariable(WLiquid, debug);
    infoFieldVariable(rhoLiquid, debug);

    // Surface tension
    const volScalarField::Internal sigma
    (
        fluid_.sigma
        (
            phaseInterface
            (
                fluid_.phases()[phaseNames().first()],
                fluid_.phases()[phaseNames().second()]
            )
        )
    );
    infoFieldVariable(sigma, debug);

    // Mole fraction of nucleating specie
    const volScalarField::Internal Xi
    (
        thermoGas.Y()[specieis().first()]*WGas/thermoGas.Wi(specieis().first())
    );

    // Saturation pressure and concentration
    const volScalarField::Internal pSat(saturationModel_->pSat(T()));
    const volScalarField::Internal cSat(pSat/p()*rhoGas/WGas);
    infoFieldVariable(pSat, debug);
    infoFieldVariable(cSat, debug);

    // Supersaturation of the nucleating specie
    const volScalarField::Internal S(Xi*p()/pSat);
    infoFieldVariable(S, true);

    // Mass, volume and diameter of one molecule in the condensed phase
    const volScalarField::Internal mMolc(WLiquid/NNA);
    const volScalarField::Internal vMolc(mMolc/rhoLiquid);
    const volScalarField::Internal dMolc(cbrt(6/pi*vMolc));
    infoFieldVariable(mMolc, debug);
    infoFieldVariable(vMolc, debug);
    infoFieldVariable(dMolc, debug);

    // Diameter of nuclei
    d_ = 4*sigma*vMolc/(k*T()*log(max(S, 1 + small)));
    infoField("d", d_);

    // ?
    const volScalarField::Internal deltaPhiStar(pi/3*sigma*sqr(d_));
    infoFieldVariable(deltaPhiStar, debug);

    // Ratio of nucleus volume to molecular volume
    const volScalarField::Internal iStar(pi/6*pow3(d_)/vMolc);
    infoFieldVariable(iStar, debug);

    // ?
    const volScalarField::Internal betaIStar1
    (
        sqrt(6*k*T()/mMolc)*sqrt((iStar + 1)/iStar)*sqr(d_/2 + dMolc/2)
    );
    infoFieldVariable(betaIStar1, debug);

    // Number-based nucleation rate; i.e., number of nuclei created per second
    // per unit volume
    const volScalarField::Internal J
    (
        betaIStar1*sqr(cSat*NNA)*exp(-deltaPhiStar/(k*T()))
       *sqrt(sigma/(k*T()))
       *2*mMolc/(pi*sqr(d_)*rhoLiquid)
    );
    infoFieldVariable(J, debug);

    // Mass transfer rate
    mDotByAlphaGas_ = J*iStar*mMolc;
    infoFieldVariable(mDotByAlphaGas_, debug);

    const volScalarField::Internal& alphaGas =
        mesh().lookupObject<volScalarField::Internal>(alphaNames().first());
    infoField("mDot", alphaGas*mDotByAlphaGas_);

    Info<< decrIndent;

    #undef infoFieldVariable
}


void Foam::fv::homogeneousCondensation::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    const label i = index(alphaNames(), eqn.psi().name());

    // !!! Note at present multiphaseEuler cannot linearise w.r.t alphaA in the
    // continuity equation for alphaB. So we can only create a linearised
    // source for this model in the gas volume-fraction equation.

    if (i == 0)
    {
        eqn -= fvm::Sp(mDotByAlphaGas_, eqn.psi());
    }
    else
    {
        massTransfer::addSup(alpha, rho, eqn);
    }
}


bool Foam::fv::homogeneousCondensation::read(const dictionary& dict)
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
