/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

inline Foam::pimpleNoLoopControl&
Foam::pimpleMultiRegionControl::pimple(const label i)
{
    return pimpleControls_[i];
}


inline Foam::solidNoLoopControl&
Foam::pimpleMultiRegionControl::solid(const label i)
{
    return solidControls_[i];
}


inline bool Foam::pimpleMultiRegionControl::pimpleTurbCorr(const label i) const
{
    return !pimpleControls_[i].turbOnFinalIterOnly() || finalPimpleIter();
}


// ************************************************************************* //
