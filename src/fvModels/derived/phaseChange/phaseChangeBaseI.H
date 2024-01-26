/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "phaseChangeBase.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::ThermoRefPair<Foam::basicThermo>&
Foam::fv::phaseChangeBase::thermos() const
{
    return thermos_;
}


inline const Foam::ThermoRefPair<Foam::fluidThermo>&
Foam::fv::phaseChangeBase::fluidThermos() const
{
    return fluidThermos_;
}


inline const Foam::ThermoRefPair<Foam::multicomponentThermo>&
Foam::fv::phaseChangeBase::specieThermos() const
{
    return specieThermos_;
}


inline const Foam::Pair<Foam::word>&
Foam::fv::phaseChangeBase::heNames() const
{
    return heNames_;
}


// ************************************************************************* //

