/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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

#include "coefficientMulticomponentMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::coefficientMulticomponentMixture<ThermoType>::
coefficientMulticomponentMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    multicomponentMixture<ThermoType>
    (
        thermoDict,
        mesh,
        phaseName
    ),
    mixture_("mixture", this->specieThermos()[0])
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const typename
Foam::coefficientMulticomponentMixture<ThermoType>::thermoMixtureType&
Foam::coefficientMulticomponentMixture<ThermoType>::cellThermoMixture
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
Foam::coefficientMulticomponentMixture<ThermoType>::thermoMixtureType&
Foam::coefficientMulticomponentMixture<ThermoType>::patchFaceThermoMixture
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


// ************************************************************************* //
