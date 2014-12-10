/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "SpecieMixture.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::SpecieMixture<MixtureType>::SpecieMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh
)
:
    MixtureType
    (
        thermoDict,
        mesh
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType>
Foam::scalar Foam::SpecieMixture<MixtureType>::nMoles
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).nMoles();
}


template<class MixtureType>
Foam::scalar Foam::SpecieMixture<MixtureType>::W
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).W();
}


template<class MixtureType>
Foam::scalar Foam::SpecieMixture<MixtureType>::Cp
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).Cp(p, T);
}


template<class MixtureType>
Foam::scalar Foam::SpecieMixture<MixtureType>::Cv
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).Cv(p, T);
}


template<class MixtureType>
Foam::scalar Foam::SpecieMixture<MixtureType>::Ha
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).Ha(p, T);
}


template<class MixtureType>
Foam::scalar Foam::SpecieMixture<MixtureType>::Hs
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).Hs(p, T);
}


template<class MixtureType>
Foam::scalar Foam::SpecieMixture<MixtureType>::Hc
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).Hc();
}


template<class MixtureType>
Foam::scalar Foam::SpecieMixture<MixtureType>::S
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).S(p, T);
}


template<class MixtureType>
Foam::scalar Foam::SpecieMixture<MixtureType>::Es
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).Es(p, T);
}


template<class MixtureType>
Foam::scalar Foam::SpecieMixture<MixtureType>::G
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).G(p, T);
}


template<class MixtureType>
Foam::scalar Foam::SpecieMixture<MixtureType>::A
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).A(p, T);
}


template<class MixtureType>
Foam::scalar Foam::SpecieMixture<MixtureType>::mu
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).mu(p, T);
}


template<class MixtureType>
Foam::scalar Foam::SpecieMixture<MixtureType>::kappa
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).kappa(p, T);
}


template<class MixtureType>
Foam::scalar Foam::SpecieMixture<MixtureType>::alphah
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).alphah(p, T);
}


template<class MixtureType>
Foam::scalar Foam::SpecieMixture<MixtureType>::rho
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->getLocalThermo(speciei).rho(p, T);
}


// ************************************************************************* //
