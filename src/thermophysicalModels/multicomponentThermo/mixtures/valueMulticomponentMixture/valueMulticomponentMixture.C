/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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

#include "valueMulticomponentMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
template<class Method, class ... Args>
Foam::scalar
Foam::valueMulticomponentMixture<ThermoType>::thermoMixtureType::
massWeighted
(
    Method psiMethod,
    const Args& ... args
) const
{
    scalar psi = 0;

    forAll(Y_, i)
    {
        psi += Y_[i]*(specieThermos_[i].*psiMethod)(args ...);
    }

    return psi;
}


template<class ThermoType>
template<class Method, class ... Args>
Foam::scalar
Foam::valueMulticomponentMixture<ThermoType>::thermoMixtureType::
harmonicMassWeighted
(
    Method psiMethod,
    const Args& ... args
) const
{
    scalar rPsi = 0;

    forAll(Y_, i)
    {
        rPsi += Y_[i]/(specieThermos_[i].*psiMethod)(args ...);
    }

    return 1/rPsi;
}


template<class ThermoType>
Foam::scalar
Foam::valueMulticomponentMixture<ThermoType>::thermoMixtureType::limit
(
    const scalar T
) const
{
    return T;
}


template<class ThermoType>
template<class Method, class ... Args>
Foam::scalar
Foam::valueMulticomponentMixture<ThermoType>::transportMixtureType::
moleWeighted
(
    Method psiMethod,
    const Args& ... args
) const
{
    scalar psi = 0;

    forAll(X_, i)
    {
        psi += X_[i]*(specieThermos_[i].*psiMethod)(args ...);
    }

    return psi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::valueMulticomponentMixture<ThermoType>::valueMulticomponentMixture
(
    const dictionary& dict
)
:
    multicomponentMixture<ThermoType>(dict),
    thermoMixture_(this->specieThermos()),
    transportMixture_(this->specieThermos())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
Foam::scalar
Foam::valueMulticomponentMixture<ThermoType>::thermoMixtureType::W() const
{
    return harmonicMassWeighted(&ThermoType::W);
}


template<class ThermoType>
Foam::scalar
Foam::valueMulticomponentMixture<ThermoType>::thermoMixtureType::rho
(
    scalar p,
    scalar T
) const
{
    return harmonicMassWeighted(&ThermoType::rho, p, T);
}


template<class ThermoType>
Foam::scalar
Foam::valueMulticomponentMixture<ThermoType>::thermoMixtureType::psi
(
    scalar p,
    scalar T
) const
{
    scalar oneByRho = 0;
    scalar psiByRho2 = 0;

    forAll(Y_, i)
    {
        const scalar rhoi = specieThermos_[i].rho(p, T);
        const scalar psii = specieThermos_[i].psi(p, T);

        oneByRho += Y_[i]/rhoi;

        if (psii > 0)
        {
            psiByRho2 += Y_[i]*psii/sqr(rhoi);
        }
    }

    return psiByRho2/sqr(oneByRho);
}


template<class ThermoType>
Foam::scalar
Foam::valueMulticomponentMixture<ThermoType>::thermoMixtureType::Hf() const
{
    return massWeighted(&ThermoType::Hf);
}


#define thermoMixtureFunction(Func)                                            \
                                                                               \
    template<class ThermoType>                                                 \
    Foam::scalar                                                               \
    Foam::valueMulticomponentMixture<ThermoType>::thermoMixtureType::Func      \
    (                                                                          \
        scalar p,                                                              \
        scalar T                                                               \
    ) const                                                                    \
    {                                                                          \
        return massWeighted(&ThermoType::Func, p, T);                          \
    }

thermoMixtureFunction(Cp)
thermoMixtureFunction(Cv)
thermoMixtureFunction(Hs)
thermoMixtureFunction(Ha)
thermoMixtureFunction(Cpv)
thermoMixtureFunction(gamma)
thermoMixtureFunction(HE)


template<class ThermoType>
Foam::scalar
Foam::valueMulticomponentMixture<ThermoType>::thermoMixtureType::THE
(
    const scalar he,
    scalar p,
    scalar T0
) const
{
    return ThermoType::T
    (
        *this,
        he,
        p,
        T0,
        &thermoMixtureType::HE,
        &thermoMixtureType::Cpv,
        &thermoMixtureType::limit
    );
}


template<class ThermoType>
Foam::scalar
Foam::valueMulticomponentMixture<ThermoType>::transportMixtureType::mu
(
    scalar p,
    scalar T
) const
{
    return moleWeighted(&ThermoType::mu, p, T);
}


template<class ThermoType>
Foam::scalar
Foam::valueMulticomponentMixture<ThermoType>::transportMixtureType::kappa
(
    scalar p,
    scalar T
) const
{
    return moleWeighted(&ThermoType::kappa, p, T);
}


template<class ThermoType>
const typename
Foam::valueMulticomponentMixture<ThermoType>::thermoMixtureType&
Foam::valueMulticomponentMixture<ThermoType>::thermoMixture
(
    const scalarFieldListSlice& Y
) const
{
    forAll(Y, i)
    {
        thermoMixture_.Y_[i] = Y[i];
    }

    return thermoMixture_;
}


template<class ThermoType>
const typename
Foam::valueMulticomponentMixture<ThermoType>::transportMixtureType&
Foam::valueMulticomponentMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice& Y
) const
{
    scalar sumX = 0;

    forAll(Y, i)
    {
        transportMixture_.X_[i] = Y[i]/this->specieThermos()[i].W();
        sumX += transportMixture_.X_[i];
    }

    forAll(Y, i)
    {
        transportMixture_.X_[i] /= sumX;
    }

    return transportMixture_;
}


template<class ThermoType>
const typename
Foam::valueMulticomponentMixture<ThermoType>::transportMixtureType&
Foam::valueMulticomponentMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice& Y,
    const thermoMixtureType&
) const
{
    return transportMixture(Y);
}


// ************************************************************************* //
