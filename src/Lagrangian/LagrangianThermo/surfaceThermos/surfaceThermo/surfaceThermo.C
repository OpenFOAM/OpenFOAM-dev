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

#include "surfaceThermo.H"
#include "LagrangianModel.H"
#include "toSubField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceThermo, 0);
    defineRunTimeSelectionTable(surfaceThermo, cloud);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::surfaceThermo::setSubMesh(const LagrangianSubMesh& subMesh) const
{
    if (subMesh.index() != subMeshIndex_)
    {
        reset();

        subMeshIndex_ = subMesh.index();
    }
}


void Foam::surfaceThermo::reset() const
{
    surface_ = false;

    resetT();
}


void Foam::surfaceThermo::implementation::resetT() const
{
    tTs_.clear();

    TmodelPtr_ = nullptr;
}


void Foam::surfaceThermo::resetDerived() const
{
    resetDerivedT();
}


void Foam::surfaceThermo::implementation::resetDerivedT() const
{
    resetDerived(tWs_);
    resetDerived(trhos_);
    resetDerived(tes_);
    resetDerived(thss_);
    resetDerived(thas_);
    resetDerived(tCps_);
    resetDerived(tCvs_);
    resetDerived(talphavs_);
    resetDerived(tkappas_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceThermo::surfaceThermo()
:
    subMeshIndex_(pTraits<uint64_t>::max),
    surface_(false)
{}


Foam::surfaceThermo::implementation::implementation
(
    const basicLagrangianSubThermo& farThermo
)
:
    farThermo_(farThermo),
    tWs_(nullptr),
    tTs_(nullptr),
    TmodelPtr_(nullptr),
    trhos_(nullptr),
    tes_(nullptr),
    thss_(nullptr),
    thas_(nullptr),
    tCps_(nullptr),
    tCvs_(nullptr),
    talphavs_(nullptr),
    tkappas_(nullptr),
    surfaceFunctionPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceThermo> Foam::surfaceThermo::New
(
    const basicLagrangianSubThermo& farThermo
)
{
    return New<surfaceThermo>(farThermo);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceThermo::~surfaceThermo()
{}


Foam::surfaceThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::IOdictionary&
Foam::surfaceThermo::implementation::properties() const
{
    return farThermo_.properties();
}


const Foam::LagrangianMesh& Foam::surfaceThermo::implementation::mesh() const
{
    return farThermo_.mesh();
}


const Foam::word&
Foam::surfaceThermo::implementation::phaseName() const
{
    return farThermo_.phaseName();
}


const Foam::surfaceFunction&
Foam::surfaceThermo::implementation::surfaceFunction() const
{
    if (!surfaceFunctionPtr_.valid())
    {
        surfaceFunctionPtr_.set
        (
            surfaceFunction::New
            (
                farThermo_.mesh(),
                farThermo_.properties()
            ).ptr()
        );
    }

    return surfaceFunctionPtr_();
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::surfaceThermo::implementation::W
(
    const LagrangianSubMesh& subMesh
) const
{
    setSubMesh(subMesh);

    if (!surface_)
    {
        return farThermo_.W(subMesh);
    }
    else
    {
        return getDerived(subMesh, tWs_);
    }
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::surfaceThermo::implementation::T
(
    const LagrangianSubMesh& subMesh
) const
{
    setSubMesh(subMesh);

    if (tTs_.valid())
    {
        return tTs_();
    }
    else
    {
        return farThermo_.T(subMesh);
    }
}


void Foam::surfaceThermo::implementation::setT
(
    const LagrangianModel& model,
    const tmp<LagrangianSubScalarSubField>& tT
)
{
    setSubMesh(tT().mesh());

    if (TmodelPtr_)
    {
        FatalErrorInFunction
            << "Lagrangian models " << TmodelPtr_->name() << " and "
            << model.name() << " both define the surface temperature"
            << exit(FatalError);
    }

    tTs_ = tT;
    tT.clear();

    TmodelPtr_ = &model;

    surface_ = true;

    resetDerived();
}


void Foam::surfaceThermo::setT
(
    const LagrangianModel& model,
    const tmp<LagrangianSubScalarField>& tT
)
{
    setT(model, toSubField(tT));
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::surfaceThermo::implementation::rho
(
    const LagrangianSubMesh& subMesh
) const
{
    setSubMesh(subMesh);

    if (!surface_)
    {
        return farThermo_.rho(subMesh);
    }
    else
    {
        return getDerived(subMesh, trhos_);
    }
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::surfaceThermo::implementation::e
(
    const LagrangianSubMesh& subMesh
) const
{
    setSubMesh(subMesh);

    if (!surface_)
    {
        return farThermo_.e(subMesh);
    }
    else
    {
        return getDerived(subMesh, tes_);
    }
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::surfaceThermo::implementation::hs
(
    const LagrangianSubMesh& subMesh
) const
{
    setSubMesh(subMesh);

    if (!surface_)
    {
        return farThermo_.hs(subMesh);
    }
    else
    {
        return getDerived(subMesh, thss_);
    }
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::surfaceThermo::implementation::ha
(
    const LagrangianSubMesh& subMesh
) const
{
    setSubMesh(subMesh);

    if (!surface_)
    {
        return farThermo_.ha(subMesh);
    }
    else
    {
        return getDerived(subMesh, thas_);
    }
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::surfaceThermo::implementation::Cp
(
    const LagrangianSubMesh& subMesh
) const
{
    setSubMesh(subMesh);

    if (!surface_)
    {
        return farThermo_.Cp(subMesh);
    }
    else
    {
        return getDerived(subMesh, tCps_);
    }
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::surfaceThermo::implementation::Cv
(
    const LagrangianSubMesh& subMesh
) const
{
    setSubMesh(subMesh);

    if (!surface_)
    {
        return farThermo_.Cv(subMesh);
    }
    else
    {
        return getDerived(subMesh, tCvs_);
    }
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::surfaceThermo::implementation::alphav
(
    const LagrangianSubMesh& subMesh
) const
{
    setSubMesh(subMesh);

    if (!surface_)
    {
        return farThermo_.alphav(subMesh);
    }
    else
    {
        return getDerived(subMesh, talphavs_);
    }
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::surfaceThermo::implementation::kappa
(
    const LagrangianSubMesh& subMesh
) const
{
    setSubMesh(subMesh);

    if (!surface_)
    {
        return farThermo_.kappa(subMesh);
    }
    else
    {
        return getDerived(subMesh, tkappas_);
    }
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::surfaceThermo::implementation::kappaNear
(
    const LagrangianSubMesh& subMesh
) const
{
    setSubMesh(subMesh);

    if (!surface_)
    {
        return farThermo_.kappa(subMesh);
    }
    else
    {
        return
            surfaceFunction().thermalConductivity
            (
                subMesh,
                farThermo_,
                static_cast<const basicLagrangianSubThermo&>(*this),
                &basicLagrangianSubThermo::kappa
            );
    }
}


// ************************************************************************* //
