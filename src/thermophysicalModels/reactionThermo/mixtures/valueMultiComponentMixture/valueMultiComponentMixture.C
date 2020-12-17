/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "valueMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::valueMultiComponentMixture<ThermoType>::valueMultiComponentMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    multiComponentMixture<ThermoType>
    (
        thermoDict,
        mesh,
        phaseName
    ),
    thermoMixture_(this->specieThermos()),
    transportMixture_(this->specieThermos())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
Foam::scalar
Foam::valueMultiComponentMixture<ThermoType>::thermoMixture::limit
(
    const scalar T
) const
{
    return T;
}


template<class ThermoType>
template<class Method, class ... Args>
Foam::scalar
Foam::valueMultiComponentMixture<ThermoType>::thermoMixture::massWeighted
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
Foam::valueMultiComponentMixture<ThermoType>::thermoMixture::
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
template<class Method, class ... Args>
Foam::scalar
Foam::valueMultiComponentMixture<ThermoType>::transportMixture::moleWeighted
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


template<class ThermoType>
Foam::scalar Foam::valueMultiComponentMixture<ThermoType>::thermoMixture::W
() const
{
    return harmonicMassWeighted(&ThermoType::W);
}


template<class ThermoType>
Foam::scalar Foam::valueMultiComponentMixture<ThermoType>::thermoMixture::rho
(
    scalar p,
    scalar T
) const
{
    return harmonicMassWeighted(&ThermoType::rho, p, T);
}


template<class ThermoType>
Foam::scalar Foam::valueMultiComponentMixture<ThermoType>::thermoMixture::psi
(
    scalar p,
    scalar T
) const
{
    scalar rho = 0;
    scalar psiByRho2 = 0;

    forAll(Y_, i)
    {
        const scalar rhoi = specieThermos_[i].rho(p, T);
        const scalar psii = specieThermos_[i].psi(p, T);

        rho += Y_[i]*rhoi;

        if (psii > 0)
        {
            psiByRho2 += Y_[i]*psii/sqr(rhoi);
        }
    }

    return sqr(rho)*psiByRho2;
}


template<class ThermoType>
Foam::scalar Foam::valueMultiComponentMixture<ThermoType>::thermoMixture::Hf
() const
{
    return massWeighted(&ThermoType::Hf);
}


#define thermoMixtureFunction(Func)                                            \
template<class ThermoType>                                                     \
Foam::scalar                                                                   \
Foam::valueMultiComponentMixture<ThermoType>::thermoMixture::Func            \
(                                                                              \
    scalar p,                                                                  \
    scalar T                                                                   \
) const                                                                        \
{                                                                              \
    return massWeighted(&ThermoType::Func, p, T);                              \
}

thermoMixtureFunction(Cp)
thermoMixtureFunction(Cv)
thermoMixtureFunction(Hs)
thermoMixtureFunction(Ha)
thermoMixtureFunction(Cpv)
thermoMixtureFunction(gamma)
thermoMixtureFunction(HE)


template<class ThermoType>
Foam::scalar Foam::valueMultiComponentMixture<ThermoType>::thermoMixture::THE
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
        &thermoMixture::HE,
        &thermoMixture::Cpv,
        &thermoMixture::limit
    );
}


template<class ThermoType>
Foam::scalar
Foam::valueMultiComponentMixture<ThermoType>::transportMixture::mu
(
    scalar p,
    scalar T
) const
{
    return moleWeighted(&ThermoType::mu, p, T);
}


template<class ThermoType>
Foam::scalar
Foam::valueMultiComponentMixture<ThermoType>::transportMixture::kappa
(
    scalar p,
    scalar T
) const
{
    return moleWeighted(&ThermoType::kappa, p, T);
}


template<class ThermoType>
const typename
Foam::valueMultiComponentMixture<ThermoType>::thermoMixtureType&
Foam::valueMultiComponentMixture<ThermoType>::cellThermoMixture
(
    const label celli
) const
{
    List<scalar>& Y = thermoMixture_.Y_;

    forAll(Y, i)
    {
        Y[i] = this->Y()[i][celli];
    }

    return thermoMixture_;
}


template<class ThermoType>
const typename
Foam::valueMultiComponentMixture<ThermoType>::thermoMixtureType&
Foam::valueMultiComponentMixture<ThermoType>::patchFaceThermoMixture
(
    const label patchi,
    const label facei
) const
{
    List<scalar>& Y = thermoMixture_.Y_;

    forAll(Y, i)
    {
        Y[i] = this->Y()[i].boundaryField()[patchi][facei];
    }

    return thermoMixture_;
}


template<class ThermoType>
const typename
Foam::valueMultiComponentMixture<ThermoType>::transportMixtureType&
Foam::valueMultiComponentMixture<ThermoType>::cellTransportMixture
(
    const label celli
) const
{
    List<scalar>& X = transportMixture_.X_;

    scalar sumX = 0;

    forAll(X, i)
    {
        X[i] = this->Y()[i][celli]/this->specieThermos()[i].W();
        sumX += X[i];
    }

    forAll(X, i)
    {
        X[i] /= sumX;
    }

    return transportMixture_;
}


template<class ThermoType>
const typename
Foam::valueMultiComponentMixture<ThermoType>::transportMixtureType&
Foam::valueMultiComponentMixture<ThermoType>::patchFaceTransportMixture
(
    const label patchi,
    const label facei
) const
{
    List<scalar>& X = transportMixture_.X_;

    scalar sumX = 0;

    forAll(X, i)
    {
        X[i] =
            this->Y()[i].boundaryField()[patchi][facei]
           /this->specieThermos()[i].W();
        sumX += X[i];
    }

    forAll(X, i)
    {
        X[i] /= sumX;
    }

    return transportMixture_;
}


// ************************************************************************* //
