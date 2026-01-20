/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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

#include "pointMasses.H"
#include "TableReader.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
    defineTypeNameAndDebug(pointMasses, 0);
    addToRunTimeSelectionTable(rigidBody, pointMasses, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::pointMasses::pointMasses
(
    const word& name,
    const dictionary& dict
)
:
    rigidBody(name, rigidBodyInertia()),
    reader_
    (
        dict.isDict("pointMasses")
      ? TableReader<point, scalar>::New
        (
            "pointMasses",
            {dimLength, dimMass},
            dict.subDict("pointMasses")
        ).ptr()
      : nullptr
    ),
    points_(),
    masses_()
{
    // Read the point masses
    const List<Tuple2<point, scalar>> pointsAndMasses =
        reader_.valid()
      ? reader_->read
        (
            {dimLength, dimMass},
            dict.subDict("pointMasses")
        )
      : TableReaders::Embedded<point, scalar>().read
        (
            {dimLength, dimMass},
            dict,
            "pointMasses"
        );
    points_.resize(pointsAndMasses.size());
    masses_.resize(pointsAndMasses.size());
    forAll(pointsAndMasses, i)
    {
        points_[i] = pointsAndMasses[i].first();
        masses_[i] = pointsAndMasses[i].second();
    }

    if (points_.size() < 4)
    {
        FatalIOErrorInFunction(dict)
            << "A " << typeName << " body requires at least four point masses"
            << exit(FatalIOError);
    }

    // Calculate the mass, the centre of mass, and the inertia tensor about the
    // centre of mass
    const scalar m = sum(masses_);
    const point c = sum(points_*masses_)/m;
    const pointField r(points_ - c);
    const symmTensor Ic = sum((symmTensor::I*(r & r) - sqr(r))*masses_);
    const scalar detIc = det(Ic);

    if (detIc < rootVSmall || detIc < rootSmall*pow3(mag(Ic)))
    {
        FatalIOErrorInFunction(dict)
            << "The inertia tensor of the " << typeName << " body '" << name
            << "' is singular. Ensure the point masses are not co-planar or"
            << " co-linear." << exit(FatalIOError);
    }

    rigidBodyInertia::operator=(rigidBodyInertia(m, c, Ic));

    Info<< *this << endl;
}


Foam::autoPtr<Foam::RBD::rigidBody> Foam::RBD::pointMasses::clone() const
{
    return autoPtr<rigidBody>(new pointMasses(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::pointMasses::~pointMasses()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::pointMasses::write(Ostream& os) const
{
    writeEntry(os, "type", type());

    // Write the point masses. Bit of a palarver to get the formatting right.
    List<Tuple2<point, scalar>> pointsAndMasses(points_.size());
    forAll(pointsAndMasses, i)
    {
        pointsAndMasses[i].first() = points_[i];
        pointsAndMasses[i].second() = masses_[i];
    }
    if (reader_.valid())
    {
        os  << indent <<  "pointMasses" << endl
            << indent << token::BEGIN_BLOCK << endl
            << incrIndent;
        reader_->write(os, {dimLength, dimMass}, pointsAndMasses);
        os  << decrIndent
            << indent << token::END_BLOCK << endl;
    }
    else
    {
        TableReaders::Embedded<point, scalar>().write
        (
            os,
            {dimLength, dimMass},
            pointsAndMasses,
            "pointMasses"
        );
    }
}


// ************************************************************************* //
