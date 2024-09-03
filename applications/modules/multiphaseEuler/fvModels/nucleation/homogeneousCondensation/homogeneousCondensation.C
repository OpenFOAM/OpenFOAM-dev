/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    singleComponentPhaseChange
    (
        name,
        modelType,
        mesh,
        dict,
        {false, false},
        {true, false}
    ),
    fluid_
    (
        mesh().lookupObject<phaseSystem>(phaseSystem::propertiesName)
    ),
    interface_
    (
        fluid_.phases()[phaseNames().first()],
        fluid_.phases()[phaseNames().second()]
    ),
    d_
    (
        IOobject
        (
            typedName(IOobject::groupName("d", interface_.name())),
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
            typedName(IOobject::groupName("mDotByAlpha", interface_.name())),
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

    const multicomponentThermo& thermoGas = specieThermos().first();

    const volScalarField& p = this->p();
    const volScalarField& T = thermoGas.T();

    const volScalarField::Internal rhoLiquid
    (
        specieThermos().valid().second()
      ? vfToVif(specieThermos().second().rhoi(specieis().second(), p, T))
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


void Foam::fv::homogeneousCondensation::correct()
{
    #define DebugField(field)                                                  \
        DebugInfo                                                              \
            << name() << ": "                                                  \
            << #field << ' ' << field.dimensions() << " min/avg/max = "        \
            << gMin(field) << '/' << gAverage(field) << '/' << gMax(field)     \
            << nl;

    using constant::mathematical::pi;
    using constant::physicoChemical::NA;
    using constant::physicoChemical::k;

    const multicomponentThermo& thermoGas = specieThermos().first();

    const volScalarField& p = this->p();
    const volScalarField& T = thermoGas.T();
    DebugField(p);
    DebugField(T);

    // Phase molecular masses and densities
    const volScalarField::Internal WGas(vfToVif(thermoGas.W()));
    const volScalarField::Internal rhoGas(vfToVif(thermoGas.rho()));
    const volScalarField::Internal WLiquid
    (
        specieThermos().valid().second()
      ? volScalarField::Internal::New
        (
            "W",
            mesh(),
            specieThermos().second().Wi(specieis().second())
        )
      : vfToVif(thermos().second().W())
    );
    const volScalarField::Internal rhoLiquid
    (
        specieThermos().valid().second()
      ? vfToVif(specieThermos().second().rhoi(specieis().second(), p, T))
      : vfToVif(thermos().second().rho())
    );
    DebugField(WGas);
    DebugField(rhoGas);
    DebugField(WLiquid);
    DebugField(rhoLiquid);

    // Surface tension
    const volScalarField::Internal sigma(interface_.sigma());
    DebugField(sigma);

    // Mole fraction of nucleating specie
    const volScalarField::Internal Xi
    (
        thermoGas.Y()[specieis().first()]*WGas/thermoGas.Wi(specieis().first())
    );

    // Saturation pressure and concentration
    const volScalarField::Internal pSat(saturationModel_->pSat(T()));
    const volScalarField::Internal cSat(pSat/p()*rhoGas/WGas);
    DebugField(pSat);
    DebugField(cSat);

    // Supersaturation of the nucleating specie
    const volScalarField::Internal S(Xi*p()/pSat);
    DebugField(S);

    // Mass, volume and diameter of one molecule in the condensed phase
    const volScalarField::Internal mMolc(WLiquid/NA);
    const volScalarField::Internal vMolc(mMolc/rhoLiquid);
    const volScalarField::Internal dMolc(cbrt(6/pi*vMolc));
    DebugField(mMolc);
    DebugField(vMolc);
    DebugField(dMolc);

    // Diameter of nuclei
    d_ = 4*sigma*vMolc/(k*T()*log(max(S, 1 + small)));
    DebugField(d_);

    // ?
    const volScalarField::Internal deltaPhiStar(pi/3*sigma*sqr(d_));
    DebugField(deltaPhiStar);

    // Ratio of nucleus volume to molecular volume
    const volScalarField::Internal iStar(pi/6*pow3(d_)/vMolc);
    DebugField(iStar);

    // ?
    const volScalarField::Internal betaIStar1
    (
        sqrt(6*k*T()/mMolc)*sqrt((iStar + 1)/iStar)*sqr(d_/2 + dMolc/2)
    );
    DebugField(betaIStar1);

    // Number-based nucleation rate; i.e., number of nuclei created per second
    // per unit volume
    const volScalarField::Internal J
    (
        betaIStar1*sqr(cSat*NA)*exp(-deltaPhiStar/(k*T()))
       *sqrt(sigma/(k*T()))
       *2*mMolc/(pi*sqr(d_)*rhoLiquid)
    );
    DebugField(J);

    // Mass transfer rate
    mDotByAlphaGas_ = J*iStar*mMolc;
    DebugField(mDotByAlphaGas_);

    #undef DebugField
}


void Foam::fv::homogeneousCondensation::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    const label i = index(alphaNames(), eqn.psi().name());

    if (i != -1)
    {
        if (i == 0)
        {
            eqn -= fvm::Sp(mDotByAlphaGas_, eqn.psi());
        }
        else
        {
            eqn += mDot() - correction(fvm::Sp(mDotByAlphaGas_, eqn.psi()));
        }
    }
    else
    {
        phaseChange::addSup(alpha, rho, eqn);
    }
}


void Foam::fv::homogeneousCondensation::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& Yi,
    fvMatrix<scalar>& eqn
) const
{
    const label i = index(phaseNames(), eqn.psi().group());

    if (i == 0 && specieis().first() != -1 && Yi.member() == specie())
    {
        eqn -= fvm::Sp(alpha*mDotByAlphaGas_, Yi);
    }
    else
    {
        singleComponentPhaseChange::addSup(alpha, rho, Yi, eqn);
    }
}


bool Foam::fv::homogeneousCondensation::read(const dictionary& dict)
{
    if (singleComponentPhaseChange::read(dict))
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
