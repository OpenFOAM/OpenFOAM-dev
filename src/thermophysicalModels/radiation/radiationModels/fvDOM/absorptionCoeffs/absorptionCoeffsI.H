/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

inline bool Foam::radiation::absorptionCoeffs::invTemp() const
{
    return  invTemp_;
}


inline Foam::scalar Foam::radiation::absorptionCoeffs::Tcommon() const
{
    return  Tcommon_;
}


inline Foam::scalar Foam::radiation::absorptionCoeffs::Tlow() const
{
    return  Tlow_;
}


inline Foam::scalar Foam::radiation::absorptionCoeffs::Thigh() const
{
    return  Thigh_;
}


inline const Foam::radiation::absorptionCoeffs::coeffArray&
Foam::radiation::absorptionCoeffs::highACoeffs() const
{
    return  highACoeffs_;
}


inline const Foam::radiation::absorptionCoeffs::coeffArray&
Foam::radiation::absorptionCoeffs::lowACoeffs() const
{
    return  lowACoeffs_;
}


// ************************************************************************* //
