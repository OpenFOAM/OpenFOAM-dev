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

#include "fluidSurfaceThermo.H"
#include "LagrangianModel.H"
#include "toSubField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidSurfaceThermo, 0);
    defineRunTimeSelectionTable(fluidSurfaceThermo, cloud);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fluidSurfaceThermo::reset() const
{
    surface_ = false;

    resetT();
    resetp();
}


void Foam::fluidSurfaceThermo::implementation::resetp() const
{
    tps_.clear();

    pModelPtr_ = nullptr;
}


void Foam::fluidSurfaceThermo::resetDerived() const
{
    resetDerivedT();
    resetDerivedp();
}


void Foam::fluidSurfaceThermo::implementation::resetDerivedp() const
{
    surfaceThermo::resetDerived(tpsis_);
    surfaceThermo::resetDerived(tmus_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidSurfaceThermo::implementation::implementation
(
    const basicLagrangianSubThermo& farThermo
)
:
    farFluidThermo_
    (
        refCast<const fluidLagrangianSubThermo>(farThermo)
    ),
    tps_(nullptr),
    pModelPtr_(nullptr),
    tpsis_(nullptr),
    tmus_(nullptr)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidSurfaceThermo> Foam::fluidSurfaceThermo::New
(
    const fluidLagrangianSubThermo& farFluidThermo
)
{
    return surfaceThermo::New<fluidSurfaceThermo>(farFluidThermo);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidSurfaceThermo::~fluidSurfaceThermo()
{}


Foam::fluidSurfaceThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::fluidSurfaceThermo::implementation::p
(
    const LagrangianSubMesh& subMesh
) const
{
    setSubMesh(subMesh);

    if (tps_.valid())
    {
        return tps_();
    }
    else
    {
        return farFluidThermo_.p(subMesh);
    }
}


void Foam::fluidSurfaceThermo::implementation::setp
(
    const LagrangianModel& model,
    const tmp<LagrangianSubScalarSubField>& tp
)
{
    setSubMesh(tp().mesh());

    if (pModelPtr_)
    {
        FatalErrorInFunction
            << "Lagrangian models " << pModelPtr_->name() << " and "
            << model.name() << " both define the surface pressure"
            << exit(FatalError);
    }

    tps_ = tp;
    tp.clear();

    pModelPtr_ = &model;

    surface_ = true;

    resetDerived();
}


void Foam::fluidSurfaceThermo::setp
(
    const LagrangianModel& model,
    const tmp<LagrangianSubScalarField>& tp
)
{
    setp(model, toSubField(tp));
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::fluidSurfaceThermo::implementation::psi
(
    const LagrangianSubMesh& subMesh
) const
{
    setSubMesh(subMesh);

    if (!surface_)
    {
        return farFluidThermo_.psi(subMesh);
    }
    else
    {
        return getDerived(subMesh, tpsis_);
    }
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::fluidSurfaceThermo::implementation::mu
(
    const LagrangianSubMesh& subMesh
) const
{
    setSubMesh(subMesh);

    if (!surface_)
    {
        return farFluidThermo_.mu(subMesh);
    }
    else
    {
        return getDerived(subMesh, tmus_);
    }
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::fluidSurfaceThermo::implementation::muNear
(
    const LagrangianSubMesh& subMesh
) const
{
    setSubMesh(subMesh);

    if (!surface_)
    {
        return farFluidThermo_.mu(subMesh);
    }
    else
    {
        return
            surfaceFunction().dynamicDiffusivity
            (
                subMesh,
                farFluidThermo_,
                static_cast<const fluidLagrangianSubThermo&>(*this),
                &fluidLagrangianSubThermo::mu
            );
    }
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::fluidSurfaceThermo::implementation::nuNear
(
    const LagrangianSubMesh& subMesh
) const
{
    if (!surface_)
    {
        return farFluidThermo_.nu(subMesh);
    }
    else
    {
        return
            surfaceFunction().kinematicDiffusivity
            (
                subMesh,
                farFluidThermo_,
                static_cast<const fluidLagrangianSubThermo&>(*this),
                &fluidLagrangianSubThermo::nu
            );
    }
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::fluidSurfaceThermo::implementation::PrNear
(
    const LagrangianSubMesh& subMesh
) const
{
    setSubMesh(subMesh);

    if (!surface_)
    {
        return farFluidThermo_.Pr(subMesh);
    }
    else
    {
        return
            muNear(subMesh)
           /surfaceFunction().dynamicDiffusivity
            (
                subMesh,
                static_cast<const basicLagrangianSubThermo&>(farFluidThermo_),
                static_cast<const basicLagrangianSubThermo&>(*this),
                &basicLagrangianSubThermo::kappaByCp
            );
    }
}


// ************************************************************************* //
