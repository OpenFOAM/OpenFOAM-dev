/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

#include "uMulticomponentMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::uMulticomponentMixture<ThermoType>::uMulticomponentMixture
(
    const dictionary& dict
)
:
    coefficientMulticomponentMixture<ThermoType>(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
Foam::scalar Foam::uMulticomponentMixture<ThermoType>::Phi
(
    const scalarFieldListSlice& Y
) const
{
    NotImplemented;
    return 1;
}


template<class ThermoType>
Foam::PtrList<Foam::volScalarField::Internal>
Foam::uMulticomponentMixture<ThermoType>::prompt
(
    const PtrList<volScalarField>& Yu
) const
{
    NotImplemented;

    PtrList<volScalarField::Internal> Yp(1);
    Yp.set(0, Yu[0]());

    return Yp;
}


// ************************************************************************* //
