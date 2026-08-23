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

#include "multicomponentSurfaceThermo.H"
#include "LagrangianModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multicomponentSurfaceThermo, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::multicomponentSurfaceThermo::reset() const
{
    surface_ = false;

    resetT();
    resetYs();
}


void Foam::multicomponentSurfaceThermo::implementation::resetYs() const
{
    Ys_.clear();

    YsModelPtr_ = nullptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multicomponentSurfaceThermo::implementation::implementation
(
    const basicLagrangianSubThermo& farThermo
)
:
    farMulticomponentThermo_
    (
        refCast<const multicomponentLagrangianSubThermo>(farThermo)
    ),
    Ys_(),
    YsModelPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multicomponentSurfaceThermo::~multicomponentSurfaceThermo()
{}


Foam::multicomponentSurfaceThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::speciesTable&
Foam::multicomponentSurfaceThermo::implementation::species() const
{
    return farMulticomponentThermo_.species();
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::multicomponentSurfaceThermo::implementation::Y
(
    const label speciei,
    const LagrangianSubMesh& subMesh
) const
{
    if (Ys_.size())
    {
        return Ys_[speciei];
    }
    else
    {
        return farMulticomponentThermo_.Y(speciei, subMesh);
    }
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::multicomponentSurfaceThermo::implementation::Y
(
    const word& specieName,
    const LagrangianSubMesh& subMesh
) const
{
    if (Ys_.size())
    {
        return Ys_[species()[specieName]];
    }
    else
    {
        return farMulticomponentThermo_.Y(specieName, subMesh);
    }
}


void Foam::multicomponentSurfaceThermo::implementation::setYs
(
    const LagrangianModel& model,
    PtrList<LagrangianSubScalarSubField>&& Ys
)
{
    if (Ys.empty()) return;

    setSubMesh(Ys[0].mesh());

    if (YsModelPtr_)
    {
        FatalErrorInFunction
            << "Lagrangian models " << YsModelPtr_->name() << " and "
            << model.name() << " both define the surface composition"
            << exit(FatalError);
    }

    Ys_.transfer(Ys);

    YsModelPtr_ = &model;

    surface_ = true;

    resetDerived();
}


void Foam::multicomponentSurfaceThermo::setYs
(
    const LagrangianModel& model,
    PtrList<LagrangianSubScalarField>&& Ys
)
{
    setYs(model, toSubField(std::move(Ys)));
}


// ************************************************************************* //
