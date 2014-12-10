/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "bladeModel.H"
#include "unitConversion.H"
#include "Tuple2.H"
#include "vector.H"
#include "IFstream.H"


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::bladeModel::readFromFile() const
{
    return fName_ != fileName::null;
}


void Foam::bladeModel::interpolateWeights
(
    const scalar& xIn,
    const List<scalar>& values,
    label& i1,
    label& i2,
    scalar& ddx
) const
{
    i2 = 0;
    label nElem = values.size();

    if (nElem == 1)
    {
        i1 = i2;
        ddx = 0.0;
        return;
    }
    else
    {
        while ((values[i2] < xIn) && (i2 < nElem))
        {
            i2++;
        }

        if (i2 == nElem)
        {
            i2 = nElem - 1;
            i1 = i2;
            ddx = 0.0;
            return;
        }
        else
        {
            i1 = i2 - 1;
            ddx = (xIn - values[i1])/(values[i2] - values[i1]);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bladeModel::bladeModel(const dictionary& dict)
:
    profileName_(),
    profileID_(),
    radius_(),
    twist_(),
    chord_(),
    fName_(fileName::null)
{
    List<Tuple2<word, vector> > data;
    if (readFromFile())
    {
        IFstream is(fName_);
        is  >> data;
    }
    else
    {
        dict.lookup("data") >> data;
    }


    if (data.size() > 0)
    {
        profileName_.setSize(data.size());
        profileID_.setSize(data.size());
        radius_.setSize(data.size());
        twist_.setSize(data.size());
        chord_.setSize(data.size());

        forAll(data, i)
        {
            profileName_[i] = data[i].first();
            profileID_[i] = -1;
            radius_[i] = data[i].second()[0];
            twist_[i] = degToRad(data[i].second()[1]);
            chord_[i] = data[i].second()[2];
        }
    }
    else
    {
        FatalErrorIn("Foam::bladeModel::bladeModel(const dictionary&)")
            << "No blade data specified" << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::bladeModel::~bladeModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::List<Foam::word>& Foam::bladeModel::profileName() const
{
    return profileName_;
}


const Foam::List<Foam::label>& Foam::bladeModel::profileID() const
{
    return profileID_;
}


const Foam::List<Foam::scalar>& Foam::bladeModel::radius() const
{
    return radius_;
}


const Foam::List<Foam::scalar>& Foam::bladeModel::twist() const
{
    return twist_;
}


const Foam::List<Foam::scalar>& Foam::bladeModel::chord() const
{
    return chord_;
}


Foam::List<Foam::label>& Foam::bladeModel::profileID()
{
    return profileID_;
}


void Foam::bladeModel::interpolate
(
    const scalar radius,
    scalar& twist,
    scalar& chord,
    label& i1,
    label& i2,
    scalar& invDr
) const
{
    interpolateWeights(radius, radius_, i1, i2, invDr);

    twist = invDr*(twist_[i2] - twist_[i1]) + twist_[i1];
    chord = invDr*(chord_[i2] - chord_[i1]) + chord_[i1];
}


// ************************************************************************* //
