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

#include "reactingMultiphaseParcelInjectionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::scalarList&
Foam::reactingMultiphaseParcelInjectionData::YGas() const
{
    return YGas_;
}


inline const Foam::scalarList&
Foam::reactingMultiphaseParcelInjectionData::YLiquid() const
{
    return YLiquid_;
}


inline const Foam::scalarList&
Foam::reactingMultiphaseParcelInjectionData::YSolid() const
{
    return YSolid_;
}


inline Foam::scalarList& Foam::reactingMultiphaseParcelInjectionData::YGas()
{
    return YGas_;
}


inline Foam::scalarList& Foam::reactingMultiphaseParcelInjectionData::YLiquid()
{
    return YLiquid_;
}


inline Foam::scalarList& Foam::reactingMultiphaseParcelInjectionData::YSolid()
{
    return YSolid_;
}


// ************************************************************************* //
