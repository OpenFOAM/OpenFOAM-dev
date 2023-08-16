/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "multicomponentThermo.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::multicomponentThermo::containsSpecie
(
    const word& specieName
) const
{
    return species().found(specieName);
}


inline void Foam::multicomponentThermo::setSpecieActive
(
    const label speciei
) const
{
    const_cast<List<bool>&>(speciesActive())[speciei] = true;
    const_cast<volScalarField&>(Y()[speciei]).writeOpt() = IOobject::AUTO_WRITE;
}


inline void Foam::multicomponentThermo::setSpecieInactive
(
    const label speciei
) const
{
    const_cast<List<bool>&>(speciesActive())[speciei] = false;
    const_cast<volScalarField&>(Y()[speciei]).writeOpt() = IOobject::NO_WRITE;
}


inline bool Foam::multicomponentThermo::solveSpecie(const label speciei) const
{
    return speciei != defaultSpecie() && speciesActive()[speciei];
}


inline Foam::volScalarField& Foam::multicomponentThermo::Y(const label speciei)
{
    return Y()[speciei];
}


inline const Foam::volScalarField& Foam::multicomponentThermo::Y
(
    const label speciei
) const
{
    return Y()[speciei];
}


inline Foam::volScalarField& Foam::multicomponentThermo::Y
(
    const word& specieName
)
{
    return Y()[species()[specieName]];
}


inline const Foam::volScalarField& Foam::multicomponentThermo::Y
(
    const word& specieName
) const
{
    return Y()[species()[specieName]];
}


inline Foam::label Foam::multicomponentThermo::specieIndex
(
    const volScalarField& Yi
) const
{
    return species()[Yi.member()];
}


inline Foam::volScalarFieldListSlicer
Foam::multicomponentThermo::implementation::Yslicer() const
{
    return volScalarFieldListSlicer(mesh(), Y_);
}


inline Foam::scalarFieldListSlice
Foam::multicomponentThermo::implementation::cellComposition
(
    const volScalarFieldListSlicer& Yslicer,
    const label celli
) const
{
    return Yslicer.slice(celli);
}


inline Foam::scalarFieldListSlice
Foam::multicomponentThermo::implementation::patchFaceComposition
(
    const volScalarFieldListSlicer& Yslicer,
    const label patchi,
    const label facei
) const
{
    return Yslicer.patchSlice(patchi, facei);
}


// ************************************************************************* //
