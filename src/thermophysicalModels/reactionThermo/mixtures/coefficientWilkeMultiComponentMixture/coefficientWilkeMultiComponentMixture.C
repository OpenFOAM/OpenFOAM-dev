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

#include "coefficientWilkeMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::coefficientWilkeMultiComponentMixture<ThermoType>::
coefficientWilkeMultiComponentMixture
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
    mixture_("mixture", this->specieThermos()[0]),
    transportMixture_(this->specieThermos())
{}


template<class ThermoType>
Foam::coefficientWilkeMultiComponentMixture<ThermoType>::transportMixture::
transportMixture
(
    const PtrList<ThermoType>& specieThermos
)
:
    specieThermos_(specieThermos),
    M_(specieThermos.size()),
    A_(specieThermos.size()),
    B_(specieThermos.size()),
    X_(specieThermos.size()),
    mu_(specieThermos.size()),
    w_(specieThermos.size()),
    muCached_(false)
{
    forAll(specieThermos_, i)
    {
        M_[i] = specieThermos[i].W();
    }

    forAll(M_, i)
    {
        forAll(M_, j)
        {
            if (i != j)
            {
                A_(i, j) = ((4/sqrt(2.0))*sqrt(1 + M_[i]/M_[j]));
                B_(i, j) = sqrt(M_[j]/M_[i]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::coefficientWilkeMultiComponentMixture<ThermoType>::transportMixture::
WilkeWeights
(
    scalar p,
    scalar T
) const
{
    forAll(mu_, i)
    {
        mu_[i] = specieThermos_[i].mu(p, T);
    }

    forAll(M_, i)
    {
        scalar sumXphi = 0;

        forAll(M_, j)
        {
            if (i != j)
            {
                const scalar phiij =
                    sqr(1 + sqrt((mu_[i]/mu_[j])*B_(i, j)))/A_(i, j);

                sumXphi += X_[j]*phiij;
            }
            else
            {
                sumXphi += X_[j];
            }
        }

        w_[i] = X_[i]/sumXphi;
    }
}


template<class ThermoType>
Foam::scalar
Foam::coefficientWilkeMultiComponentMixture<ThermoType>::transportMixture::
mu
(
    scalar p,
    scalar T
) const
{
    WilkeWeights(p, T);

    scalar mu = 0;
    forAll(w_, i)
    {
        mu += w_[i]*mu_[i];
    }

    return mu;
}


template<class ThermoType>
Foam::scalar
Foam::coefficientWilkeMultiComponentMixture<ThermoType>::transportMixture::
kappa
(
    scalar p,
    scalar T
) const
{
    if (!muCached_)
    {
        WilkeWeights(p, T);
    }

    scalar kappa = 0;
    forAll(w_, i)
    {
        kappa += w_[i]*specieThermos_[i].kappa(p, T);
    }

    return kappa;
}


template<class ThermoType>
const typename
Foam::coefficientWilkeMultiComponentMixture<ThermoType>::thermoMixtureType&
Foam::coefficientWilkeMultiComponentMixture<ThermoType>::cellThermoMixture
(
    const label celli
) const
{
    mixture_ = this->Y()[0][celli]*this->specieThermos()[0];

    for (label i=1; i<this->Y().size(); i++)
    {
        mixture_ += this->Y()[i][celli]*this->specieThermos()[i];
    }

    return mixture_;
}


template<class ThermoType>
const typename
Foam::coefficientWilkeMultiComponentMixture<ThermoType>::thermoMixtureType&
Foam::coefficientWilkeMultiComponentMixture<ThermoType>::patchFaceThermoMixture
(
    const label patchi,
    const label facei
) const
{
    mixture_ =
        this->Y()[0].boundaryField()[patchi][facei]
       *this->specieThermos()[0];

    for (label i=1; i<this->Y().size(); i++)
    {
        mixture_ +=
            this->Y()[i].boundaryField()[patchi][facei]
           *this->specieThermos()[i];
    }

    return mixture_;
}


template<class ThermoType>
const typename
Foam::coefficientWilkeMultiComponentMixture<ThermoType>::transportMixtureType&
Foam::coefficientWilkeMultiComponentMixture<ThermoType>::cellTransportMixture
(
    const label celli
) const
{
    transportMixture_.muCached_ = false;

    scalarList& X = transportMixture_.X_;

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
Foam::coefficientWilkeMultiComponentMixture<ThermoType>::transportMixtureType&
Foam::coefficientWilkeMultiComponentMixture<ThermoType>::
patchFaceTransportMixture
(
    const label patchi,
    const label facei
) const
{
    transportMixture_.muCached_ = false;

    scalarList& X = transportMixture_.X_;

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

template<class ThermoType>
const typename
Foam::coefficientWilkeMultiComponentMixture<ThermoType>::transportMixtureType&
Foam::coefficientWilkeMultiComponentMixture<ThermoType>::
cellTransportMixture
(
    const label celli,
    const thermoMixtureType& thermoMixture
) const
{
    cellTransportMixture(celli);
    transportMixture_.muCached_ = true;
    return transportMixture_;
}


template<class ThermoType>
const typename
Foam::coefficientWilkeMultiComponentMixture<ThermoType>::transportMixtureType&
Foam::coefficientWilkeMultiComponentMixture<ThermoType>::
patchFaceTransportMixture
(
    const label patchi,
    const label facei,
    const thermoMixtureType& thermoMixture
) const
{
    patchFaceTransportMixture(patchi, facei);
    transportMixture_.muCached_ = true;
    return transportMixture_;
}


// ************************************************************************* //
