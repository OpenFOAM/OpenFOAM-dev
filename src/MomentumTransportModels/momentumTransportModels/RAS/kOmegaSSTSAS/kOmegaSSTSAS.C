/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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

#include "kOmegaSSTSAS.H"
#include "fvcLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> kOmegaSSTSAS<BasicMomentumTransportModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    volScalarField::Internal L
    (
        sqrt(this->k_())/(pow025(this->betaStar_)*this->omega_())
    );

    volScalarField::Internal Lvk
    (
        max
        (
            kappa_*sqrt(S2)
           /(
                mag(fvc::laplacian(this->U_))()()
              + dimensionedScalar
                (
                    "rootVSmall",
                    dimensionSet(0, -1, -1, 0, 0),
                    rootVSmall
                )
            ),
            Cs_*sqrt(kappa_*zeta2_/(beta/this->betaStar_ - gamma))*delta()()
        )
    );

    return fvm::Su
    (
        this->alpha_()*this->rho_()
       *min
        (
            max
            (
                zeta2_*kappa_*S2*sqr(L/Lvk)
              - (2*C_/sigmaPhi_)*this->k_()
               *max
                (
                    magSqr(fvc::grad(this->omega_)()())/sqr(this->omega_()),
                    magSqr(fvc::grad(this->k_)()())/sqr(this->k_())
                ),
                dimensionedScalar(dimensionSet(0, 0, -2, 0, 0), 0)
            ),
            // Limit SAS production of omega for numerical stability,
            // particularly during start-up
            this->omega_()/(0.1*this->omega_.time().deltaT())
        ),
        this->omega_
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
kOmegaSSTSAS<BasicMomentumTransportModel>::kOmegaSSTSAS
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
    kOmegaSST<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.11
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    zeta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "zeta2",
            this->coeffDict_,
            3.51
        )
    ),
    sigmaPhi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaPhi",
            this->coeffDict_,
            2.0/3.0
        )
    ),
    C_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C",
            this->coeffDict_,
            2
        )
    ),

    delta_
    (
        LESdelta::New
        (
            IOobject::groupName("delta", alphaRhoPhi.group()),
            *this,
            this->coeffDict_
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool kOmegaSSTSAS<BasicMomentumTransportModel>::read()
{
    if (kOmegaSST<BasicMomentumTransportModel>::read())
    {
        Cs_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());
        sigmaPhi_.readIfPresent(this->coeffDict());
        zeta2_.readIfPresent(this->coeffDict());
        C_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
