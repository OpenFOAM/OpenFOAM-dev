/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "thermoSurfaceFilm.H"
#include "filmRadiationModel.H"
#include "heatTransferModel.H"
#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::fluidThermo&
Foam::thermoSurfaceFilm::primaryThermo() const
{
    return primaryThermo_;
}


inline const Foam::volScalarField::Internal&
Foam::thermoSurfaceFilm::hSp() const
{
    return hSp_;
}


inline const Foam::volScalarField&
Foam::thermoSurfaceFilm::hSpPrimary() const
{
    return hSpPrimary_;
}


inline const Foam::volScalarField&
Foam::thermoSurfaceFilm::TPrimary() const
{
    return TPrimary_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::thermoSurfaceFilm::YPrimary() const
{
    return YPrimary_;
}


inline const Foam::surfaceFilmModels::heatTransferModel&
Foam::thermoSurfaceFilm::htcs() const
{
    return htcs_();
}


inline const Foam::surfaceFilmModels::heatTransferModel&
Foam::thermoSurfaceFilm::htcw() const
{
    return htcw_();
}


inline const Foam::surfaceFilmModels::phaseChangeModel&
Foam::thermoSurfaceFilm::phaseChange() const
{
    return phaseChange_();
}


inline const Foam::surfaceFilmModels::radiationModel&
Foam::thermoSurfaceFilm::radiation() const
{
    return radiation_();
}


// ************************************************************************* //
