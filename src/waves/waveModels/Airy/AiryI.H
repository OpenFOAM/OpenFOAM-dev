/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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

#include "Airy.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline Foam::waveModels::AiryCoeffs Foam::waveModels::Airy::coeffs
(
    const scalar t
) const
{
    return AiryCoeffs(depth_, amplitude_->value(t), length_, g());
}


inline Foam::waveModels::AiryCoeffs Foam::waveModels::Airy::coeffs() const
{
    return AiryCoeffs(depth_, amplitude(), length_, g());
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar Foam::waveModels::Airy::depth() const
{
    return depth_;
}


inline Foam::scalar Foam::waveModels::Airy::amplitude(const scalar t) const
{
    return amplitude_->value(t);
}


inline Foam::scalar Foam::waveModels::Airy::amplitude() const
{
    return amplitude_->value(great);
}


inline Foam::scalar Foam::waveModels::Airy::length() const
{
    return length_;
}


inline Foam::scalar Foam::waveModels::Airy::phase() const
{
    return phase_;
}


// ************************************************************************* //
