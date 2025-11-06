/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "GaussianbXiIgnition.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(GaussianbXiIgnition, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            GaussianbXiIgnition,
            dictionary
        );
    }
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::GaussianbXiIgnition::calcVk()
{
    const labelList& cells = zone_.zone();
    const scalarField& V(mesh().V());
    const vectorField& C(mesh().C());

    const scalar diameter = diameter_.value();
    const vector& position = position_.value();

    Vk_.setSize(cells.size());

    forAll(cells, i)
    {
        const label celli = cells[i];
        Vk_[i] = V[celli]*exp(-magSqr(C[celli] - position)/sqr(diameter));
    }

    const scalar sumVk = gSum(Vk_);
    const scalar Vg =
        kernelShape_->Dcorr().value()
       *pow(sqrt(pi)*diameter, mesh().nGeometricD());

    forAll(Vk_, i)
    {
        Vk_[i] *= Vg/sumVk;
    }
}


void Foam::fv::GaussianbXiIgnition::readCoeffs(const dictionary& dict)
{
    position_.read(dict);
    diameter_.read(dict);

    rate_.reset
    (
        Function1<scalar>::New
        (
            "rate",
            mesh().time().userUnits(),
            dimRate,
            dict
        ).ptr()
    );

    calcVk();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::GaussianbXiIgnition::GaussianbXiIgnition
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    bXiTimedIgnition(name, modelType, mesh, dict),
    zone_(mesh, coeffs(dict)),
    kernelShape_(kernelShape::New(mesh, coeffs(dict))),
    position_("position", dimLength, Zero),
    diameter_("diameter", dimLength, 0)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::GaussianbXiIgnition::addSup
(
    const volScalarField& rho,
    const volScalarField& b,
    fvMatrix<scalar>& eqn
) const
{
    if (!igniting()) return;

    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    scalarField& Sp = eqn.diag();
    const labelList& cells = zone_.zone();
    const scalar rate = rate_->value(ignRelTime(mesh().time().value()));

    forAll(cells, i)
    {
        const label celli = cells[i];
        Sp[celli] -= Vk_[i]*rho[celli]*rate;
    }
}


void Foam::fv::GaussianbXiIgnition::topoChange
(
    const polyTopoChangeMap& map
)
{
    zone_.topoChange(map);
    calcVk();
}


void Foam::fv::GaussianbXiIgnition::mapMesh(const polyMeshMap& map)
{
    zone_.mapMesh(map);
    calcVk();
}


void Foam::fv::GaussianbXiIgnition::distribute
(
    const polyDistributionMap& map
)
{
    zone_.distribute(map);
    calcVk();
}


bool Foam::fv::GaussianbXiIgnition::movePoints()
{
    zone_.movePoints();
    calcVk();
    return true;
}


bool Foam::fv::GaussianbXiIgnition::read(const dictionary& dict)
{
    if (bXiIgnition::read(dict))
    {
        zone_.read(coeffs(dict));
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }

    return false;
}


// ************************************************************************* //
