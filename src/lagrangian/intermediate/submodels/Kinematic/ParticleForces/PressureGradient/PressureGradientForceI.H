/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::interpolation<Foam::vector>&
Foam::PressureGradientForce<CloudType>::DUcDtInterp() const
{
    if (!DUcDtInterpPtr_.valid())
    {
        FatalErrorIn
        (
            "inline const Foam::interpolation<Foam::vector>&"
            "Foam::PressureGradientForce<CloudType>::DUcDtInterp() const"
        )   << "Carrier phase DUcDt interpolation object not set"
            << abort(FatalError);
    }

    return DUcDtInterpPtr_();
}


// ************************************************************************* //
