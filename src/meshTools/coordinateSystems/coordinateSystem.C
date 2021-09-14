/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "axesRotation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coordinateSystem, 0);
    defineRunTimeSelectionTable(coordinateSystem, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const point& origin
)
:
    name_(name),
    origin_(origin),
    R_(new axesRotation(sphericalTensor::I))
{}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const point& origin,
    const coordinateRotation& cr
)
:
    name_(name),
    origin_(origin),
    R_(cr.clone())
{}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const point& origin,
    const vector& axis,
    const vector& dirn
)
:
    name_(name),
    origin_(origin),
    R_(new axesRotation(axis, dirn))
{}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    origin_(dict.lookup("origin")),
    R_(coordinateRotation::New(dict.subDict("coordinateRotation")).ptr())
{}


Foam::coordinateSystem::coordinateSystem(const coordinateSystem& cs)
:
    name_(cs.name_),
    origin_(cs.origin_),
    R_(cs.R_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordinateSystem::~coordinateSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::coordinateSystem::localToGlobal
(
    const vector& local,
    bool translate
) const
{
    if (translate)
    {
        return (R_->transform(local)) + origin_;
    }
    else
    {
        return R_->transform(local);
    }
}


Foam::tmp<Foam::vectorField> Foam::coordinateSystem::localToGlobal
(
    const vectorField& local,
    bool translate
) const
{
    if (translate)
    {
        return (R_->transform(local)) + origin_;
    }
    else
    {
        return R_->transform(local);
    }
}


Foam::vector Foam::coordinateSystem::globalToLocal
(
    const vector& global,
    bool translate
) const
{
    if (translate)
    {
        return R_->invTransform(global - origin_);
    }
    else
    {
        return R_->invTransform(global);
    }
}


Foam::tmp<Foam::vectorField> Foam::coordinateSystem::globalToLocal
(
    const vectorField& global,
    bool translate
) const
{
    if (translate)
    {
        return R_->invTransform(global - origin_);
    }
    else
    {
        return R_->invTransform(global);
    }
}


void Foam::coordinateSystem::write(Ostream& os) const
{
    os  << type() << " origin: " << origin() << nl;
    R_->write(os);
}


void Foam::coordinateSystem::writeDict(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os  << indent << name_ << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    writeEntry(os, "type", type());
    writeEntry(os, "origin", origin_);

    os  << indent << "coordinateRotation" << nl
        << indent << token::BEGIN_BLOCK << incrIndent << nl;
    writeEntry(os, "type", R_->type());
    R_->write(os);
    os  << decrIndent << indent << token::END_BLOCK << endl;

    if (subDict)
    {
        os  << decrIndent << indent << token::END_BLOCK << endl;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::coordinateSystem::operator=(const coordinateSystem& cs)
{
    name_ = cs.name_;
    origin_ = cs.origin_;
    R_ = cs.R_->clone();
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const coordinateSystem& cs)
{
    cs.writeDict(os);
    os.check("Ostream& operator<<(Ostream&, const coordinateSystem&");
    return os;
}


// ************************************************************************* //
