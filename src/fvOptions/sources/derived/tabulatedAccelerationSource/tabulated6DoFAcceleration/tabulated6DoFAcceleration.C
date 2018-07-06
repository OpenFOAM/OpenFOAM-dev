/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "tabulated6DoFAcceleration.H"
#include "Tuple2.H"
#include "IFstream.H"
#include "interpolateSplineXY.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tabulated6DoFAcceleration, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tabulated6DoFAcceleration::tabulated6DoFAcceleration
(
    const dictionary& accelerationCoeffs,
    const Time& runTime
)
:
    time_(runTime),
    accelerationCoeffs_(accelerationCoeffs)
{
    read(accelerationCoeffs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tabulated6DoFAcceleration::~tabulated6DoFAcceleration()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tabulated6DoFAcceleration::accelerationVectors
Foam::tabulated6DoFAcceleration::acceleration() const
{
    scalar t = time_.value();

    if (t < times_[0])
    {
        FatalErrorInFunction
            << "current time (" << t
            << ") is less than the minimum in the data table ("
            << times_[0] << ')'
            << exit(FatalError);
    }

    if (t > times_.last())
    {
        FatalErrorInFunction
            << "current time (" << t
            << ") is greater than the maximum in the data table ("
            << times_.last() << ')'
            << exit(FatalError);
    }

    accelerationVectors avs = interpolateSplineXY
    (
        t,
        times_,
        values_
    );

    InfoInFunction
        << "Time = " << t << " accelerations: " << avs << endl;

    return avs;
}


bool Foam::tabulated6DoFAcceleration::read
(
    const dictionary& accelerationCoeffs
)
{
    accelerationCoeffs_ = accelerationCoeffs;

    // If the timeDataFileName has changed read the file

    fileName newTimeDataFileName
    (
        fileName(accelerationCoeffs_.lookup("timeDataFileName")).expand()
    );

    if (newTimeDataFileName != timeDataFileName_)
    {
        timeDataFileName_ = newTimeDataFileName;

        IFstream dataStream(timeDataFileName_);

        if (dataStream.good())
        {
            List<Tuple2<scalar, accelerationVectors>> timeValues
            (
                dataStream
            );

            times_.setSize(timeValues.size());
            values_.setSize(timeValues.size());

            forAll(timeValues, i)
            {
                times_[i] = timeValues[i].first();
                values_[i] = timeValues[i].second();
            }
        }
        else
        {
            FatalErrorInFunction
                << "Cannot open time data file " << timeDataFileName_
                << exit(FatalError);
        }
    }

    return true;
}


// ************************************************************************* //
