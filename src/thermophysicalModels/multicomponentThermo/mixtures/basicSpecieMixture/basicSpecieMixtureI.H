/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

inline const Foam::speciesTable& Foam::basicSpecieMixture::species() const
{
    return species_;
}


inline bool Foam::basicSpecieMixture::contains(const word& specieName) const
{
    return species_.found(specieName);
}


inline Foam::label Foam::basicSpecieMixture::defaultSpecie() const
{
    return defaultSpecieIndex_;
}


inline bool Foam::basicSpecieMixture::active(label speciei) const
{
    return active_[speciei];
}


inline const Foam::List<bool>& Foam::basicSpecieMixture::active() const
{
    return active_;
}


inline void Foam::basicSpecieMixture::setActive(label speciei) const
{
    active_[speciei] = true;
    const_cast<volScalarField&>(Y_[speciei]).writeOpt() = IOobject::AUTO_WRITE;
}


inline void Foam::basicSpecieMixture::setInactive(label speciei) const
{
    active_[speciei] = false;
    const_cast<volScalarField&>(Y_[speciei]).writeOpt() = IOobject::NO_WRITE;
}


inline bool Foam::basicSpecieMixture::solve(label speciei) const
{
    return speciei != defaultSpecieIndex_ && active_[speciei];
}


inline Foam::PtrList<Foam::volScalarField>& Foam::basicSpecieMixture::Y()
{
    return Y_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basicSpecieMixture::Y() const
{
    return Y_;
}


inline Foam::volScalarField& Foam::basicSpecieMixture::Y(const label i)
{
    return Y_[i];
}


inline const Foam::volScalarField& Foam::basicSpecieMixture::Y
(
    const label i
) const
{
    return Y_[i];
}


inline Foam::volScalarField& Foam::basicSpecieMixture::Y
(
    const word& specieName
)
{
    return Y_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basicSpecieMixture::Y
(
    const word& specieName
) const
{
    return Y_[species_[specieName]];
}


inline Foam::label Foam::basicSpecieMixture::index
(
    const volScalarField& Yi
) const
{
    return species_[Yi.member()];
}


// ************************************************************************* //
