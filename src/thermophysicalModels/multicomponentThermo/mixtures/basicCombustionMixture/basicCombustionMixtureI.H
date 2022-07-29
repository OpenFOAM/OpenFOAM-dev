/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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

inline bool Foam::basicCombustionMixture::contains
(
    const word& specieName
) const
{
    return species_.found(specieName);
}


inline Foam::volScalarField& Foam::basicCombustionMixture::Y
(
    const word& specieName
)
{
    return Y_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basicCombustionMixture::Y
(
    const word& specieName
) const
{
    return Y_[species_[specieName]];
}


inline Foam::scalar Foam::basicCombustionMixture::fres
(
    const scalar ft,
    const scalar stoicRatio
) const
{
    return max(ft - (scalar(1) - ft)/stoicRatio, scalar(0));
}


inline Foam::tmp<Foam::volScalarField> Foam::basicCombustionMixture::fres
(
    const volScalarField& ft,
    const dimensionedScalar& stoicRatio
) const
{
    return max(ft - (scalar(1) - ft)/stoicRatio.value(), scalar(0));
}


// ************************************************************************* //
