/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "FluidSurfaceThermo.H"
#include "LagrangianSubFields.H"
#include "toSubField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BaseThermo>
void Foam::FluidSurfaceThermo<BaseThermo>::correct
(
    const LagrangianSubMesh& subMesh
) const
{
    using namespace dimensions;

    if (!this->surface_) return;

    auto doSubField = [&](const tmp<LagrangianSubScalarSubField>& tf)
    {
        return tf.valid() && isNull(tf());
    };
    auto doField = [&](const tmp<LagrangianSubScalarField>& tf)
    {
        return tf.valid() && isNull(tf());
    };
    const bool doW = doField(this->tWs_);
    const bool doRho = doSubField(this->trhos_);
    const bool doE = doSubField(this->tes_);
    const bool doHs = doField(this->thss_);
    const bool doHa = doField(this->thas_);
    const bool doCp = doField(this->tCps_);
    const bool doCv = doSubField(this->tCvs_);
    const bool doAlphav = doField(this->talphavs_);
    const bool doPsi = doSubField(this->tpsis_);
    const bool doKappa = doSubField(this->tkappas_);
    const bool doMu = doSubField(this->tmus_);

    const bool doThermo =
        doW || doRho || doE || doHs || doHa
     || doCp || doCv || doPsi || doAlphav;
    const bool doTransport = doKappa || doMu;

    if (!doThermo && !doTransport) return;

    auto initSubField = [&](const word& name, const dimensionSet& dims)
    {
        return toSubField<scalar, LagrangianSubMesh>(name, subMesh, dims);
    };
    auto initField = [&](const word& name, const dimensionSet& dims)
    {
        return LagrangianSubScalarField::New(name, subMesh, dims);
    };
    if (doW) this->tWs_ = initField("W", mass/moles);
    if (doRho) this->trhos_ = initSubField("rho", density);
    if (doE) this->tes_ = initSubField("e", specificEnergy);
    if (doHs) this->thss_ = initField("hs", specificEnergy);
    if (doHa) this->thas_ = initField("ha", specificEnergy);
    if (doCp) this->tCps_ = initField("Cp", specificHeatCapacity);
    if (doCv) this->tCvs_ = initSubField("Cv", specificHeatCapacity);
    if (doPsi) this->tpsis_ = initSubField("psi", compressibility);
    if (doAlphav) this->talphavs_ = initField("alphav", inv(temperature));
    if (doKappa) this->tkappas_ = initSubField("kappa", thermalConductivity);
    if (doMu) this->tmus_ = initSubField("mu", dynamicViscosity);

    tmp<LagrangianSubScalarSubField> tps = this->p(subMesh);
    const LagrangianSubScalarSubField& ps = tps();
    tmp<LagrangianSubScalarSubField> tTs = this->T(subMesh);
    const LagrangianSubScalarSubField& Ts = tTs();

    auto Yslicer = this->Yslicer(subMesh);

    forAll(subMesh, subi)
    {
        const scalar p = ps[subi], T = Ts[subi];

        auto composition = this->elementComposition(Yslicer, subi);

        const typename BaseThermo::mixtureType::thermoMixtureType&
            thermoMixture = this->thermoMixture(composition);

        if (doW) this->tWs_.ref()[subi] = thermoMixture.W();
        if (doRho) this->trhos_.ref()[subi] = thermoMixture.rho(p, T);
        if (doE) this->tes_.ref()[subi] = thermoMixture.es(p, T);
        if (doHs) this->thss_.ref()[subi] = thermoMixture.hs(p, T);
        if (doHa) this->thas_.ref()[subi] = thermoMixture.ha(p, T);
        if (doCp) this->tCps_.ref()[subi] = thermoMixture.Cp(p, T);
        if (doCv) this->tCvs_.ref()[subi] = thermoMixture.Cv(p, T);
        if (doPsi) this->tpsis_.ref()[subi] = thermoMixture.psi(p, T);
        if (doAlphav) this->talphavs_.ref()[subi] = thermoMixture.alphav(p, T);

        if (!doTransport) continue;

        const typename BaseThermo::mixtureType::transportMixtureType&
            transportMixture =
            this->transportMixture(composition, thermoMixture);

        if (doKappa) this->tkappas_.ref()[subi] = transportMixture.kappa(p, T);
        if (doMu) this->tmus_.ref()[subi] = transportMixture.mu(p, T);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::FluidSurfaceThermo<BaseThermo>::FluidSurfaceThermo
(
    const basicLagrangianSubThermo& farThermo
)
:
    BaseThermo(farThermo)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::FluidSurfaceThermo<BaseThermo>::~FluidSurfaceThermo()
{}


// ************************************************************************* //
