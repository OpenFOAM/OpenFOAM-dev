/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

#include "kOmegaSSTDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
kOmegaSSTDES<BasicMomentumTransportModel>::Lt() const
{
    return sqrt(this->k_())/(this->betaStar_*this->omega_());
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal> kOmegaSSTDES<BasicMomentumTransportModel>::FDES
(
    const volScalarField::Internal& F1,
    const volScalarField::Internal& F2
) const
{
    switch (FSST_)
    {
        case 0:
            return max(Lt()/(CDES_*this->delta()()), scalar(1));
        case 1:
            return max(Lt()*(1 - F1)/(CDES_*this->delta()()), scalar(1));
        case 2:
            return max(Lt()*(1 - F2)/(CDES_*this->delta()()), scalar(1));
        default:
            FatalErrorInFunction
                << "Incorrect FSST = " << FSST_ << ", should be 0, 1 or 2"
                << exit(FatalError);
            return F1;
    }
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
kOmegaSSTDES<BasicMomentumTransportModel>::epsilonByk
(
    const volScalarField::Internal& F1,
    const volScalarField::Internal& F2
) const
{
    return this->betaStar_*this->omega_()*FDES(F1, F2);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
kOmegaSSTDES<BasicMomentumTransportModel>::kOmegaSSTDES
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    kOmegaSST
    <
        LESeddyViscosity<BasicMomentumTransportModel>,
        BasicMomentumTransportModel
    >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    CDES_("CDES", this->coeffDict(), 0.61),
    FSST_(this->coeffDict().lookupOrDefault("FSST", 2))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool kOmegaSSTDES<BasicMomentumTransportModel>::read()
{
    if
    (
        kOmegaSST
        <
            LESeddyViscosity<BasicMomentumTransportModel>,
            BasicMomentumTransportModel
        >::read()
    )
    {
        CDES_.readIfPresent(this->coeffDict());
        this->coeffDict().readIfPresent("FSST", FSST_);

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
