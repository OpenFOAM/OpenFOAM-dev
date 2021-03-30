/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

#include "SmagorinskyZhang.H"
#include "fvModels.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
SmagorinskyZhang<BasicMomentumTransportModel>::SmagorinskyZhang
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& type
)
:
    Smagorinsky<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        type
    ),

    gasTurbulencePtr_(nullptr),

    Cmub_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmub",
            this->coeffDict_,
            0.6
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool SmagorinskyZhang<BasicMomentumTransportModel>::read()
{
    if (Smagorinsky<BasicMomentumTransportModel>::read())
    {
        Cmub_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
const PhaseCompressibleMomentumTransportModel
<
    typename BasicMomentumTransportModel::transportModel
>&
SmagorinskyZhang<BasicMomentumTransportModel>::gasTurbulence() const
{
    if (!gasTurbulencePtr_)
    {
        const volVectorField& U = this->U_;

        const phaseModel& liquid = refCast<const phaseModel>(this->transport());
        const phaseSystem& fluid = liquid.fluid();
        const phaseModel& gas = fluid.otherPhase(liquid);

        gasTurbulencePtr_ =
           &U.db().lookupObject
            <
               PhaseCompressibleMomentumTransportModel<transportModel>
            >
            (
                IOobject::groupName
                (
                    momentumTransportModel::typeName,
                    gas.name()
                )
            );
    }

    return *gasTurbulencePtr_;
}


template<class BasicMomentumTransportModel>
void SmagorinskyZhang<BasicMomentumTransportModel>::correctNut()
{
    const PhaseCompressibleMomentumTransportModel<transportModel>&
        gasTurbulence = this->gasTurbulence();

    const phaseModel& liquid = refCast<const phaseModel>(this->transport());
    const phaseSystem& fluid = liquid.fluid();
    const phaseModel& gas = fluid.otherPhase(liquid);

    volScalarField k(this->k(fvc::grad(this->U_)));

    this->nut_ =
        this->Ck_*sqrt(k)*this->delta()
      + Cmub_*gas.d()*gasTurbulence.alpha()
       *(mag(this->U_ - gasTurbulence.U()));

    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
