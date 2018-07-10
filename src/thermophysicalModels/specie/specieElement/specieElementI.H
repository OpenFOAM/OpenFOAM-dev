/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "token.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::specieElement::specieElement()
{}


inline Foam::specieElement::specieElement(const word& name, const label nAtoms)
:
    name_(name),
    nAtoms_(nAtoms)
{}


inline Foam::specieElement::specieElement(Istream& is)
:
    name_(is),
    nAtoms_(readLabel(is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::specieElement::name() const
{
    return name_;
}


inline Foam::word& Foam::specieElement::name()
{
    return name_;
}


inline Foam::label Foam::specieElement::nAtoms() const
{
    return nAtoms_;
}


inline Foam::label& Foam::specieElement::nAtoms()
{
    return nAtoms_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::specieElement::operator==(const specieElement& se) const
{
    return
    (
        nAtoms_ == se.nAtoms_
     && name_ == se.name_
    );
}


inline bool Foam::specieElement::operator!=(const specieElement& se) const
{
    return !operator==(se);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<(Ostream& os, const specieElement& se)
{
    os  << se.name() << token::SPACE << se.nAtoms();
    return os;
}


// ************************************************************************* //
