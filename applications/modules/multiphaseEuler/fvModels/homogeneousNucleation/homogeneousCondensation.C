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
    saturationModel_.reset
    (
        saturationPressureModel::New("pSat", dict).ptr()
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::fv::homogeneousCondensation::dAndMDotByAlphaSolution() const
{
    #define infoFieldVariable(field, print) infoField(#field, field, print)

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
    const volScalarField::Internal sigma(this->sigma());
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

    // Mass and diameter of one molecule in the liquid
    const volScalarField::Internal mMolc(WLiquid/NNA);
    const volScalarField::Internal dMolc(cbrt(6/pi*(mMolc/rhoLiquid)));
    infoFieldVariable(mMolc, debug);
    infoFieldVariable(dMolc, debug);

    // Diameter and mass of a nucleus
    tmp<volScalarField::Internal> td =
        4*sigma*mMolc/rhoLiquid/(k*T()*log(max(S, 1 + small)));
    const volScalarField::Internal d = td();
    const volScalarField::Internal m(pi/6*pow3(d)*rhoLiquid);
    infoFieldVariable(d, true);
    infoField("m/mMolc", m/mMolc, debug);

    // Free energy cost of a nucleus
    const volScalarField::Internal deltaPhiStar(pi/3*sigma*sqr(d));
    infoFieldVariable(deltaPhiStar, debug);

    // ?
    const volScalarField::Internal betaIStar1
    (
        sqrt(6*k*T()*(1/mMolc + 1/m))*sqr(dMolc/2 + d/2)
    );
    infoFieldVariable(betaIStar1, debug);

    // Number-based nucleation rate; i.e., number of nuclei created per second
    // per unit volume
    const volScalarField::Internal J
    (
        betaIStar1
       *sqr(cSat*NNA)
       *exp(-deltaPhiStar/(k*T()))
       *sqrt(sigma/(k*T()))
       *2*mMolc/(pi*sqr(d)*rhoLiquid)
    );
    infoFieldVariable(J, debug);

    return Pair<tmp<volScalarField::Internal>>(td, J*m);

    #undef infoFieldVariable
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
    homogeneousNucleation(name, modelType, mesh, dict),
    saturationModel_(nullptr)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::homogeneousCondensation::read(const dictionary& dict)
{
    if (homogeneousNucleation::read(dict))
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
