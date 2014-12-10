/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "tabulated6DoFMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "Tuple2.H"
#include "IFstream.H"
#include "interpolateSplineXY.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(tabulated6DoFMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        tabulated6DoFMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::tabulated6DoFMotion::tabulated6DoFMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime)
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::tabulated6DoFMotion::~tabulated6DoFMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::tabulated6DoFMotion::transformation() const
{
    scalar t = time_.value();

    if (t < times_[0])
    {
        FatalErrorIn
        (
            "solidBodyMotionFunctions::tabulated6DoFMotion::transformation()"
        )   << "current time (" << t
            << ") is less than the minimum in the data table ("
            << times_[0] << ')'
            << exit(FatalError);
    }

    if (t > times_.last())
    {
        FatalErrorIn
        (
            "solidBodyMotionFunctions::tabulated6DoFMotion::transformation()"
        )   << "current time (" << t
            << ") is greater than the maximum in the data table ("
            << times_.last() << ')'
            << exit(FatalError);
    }

    translationRotationVectors TRV = interpolateSplineXY
    (
        t,
        times_,
        values_
    );

    // Convert the rotational motion from deg to rad
    TRV[1] *= pi/180.0;

    quaternion R(TRV[1].x(), TRV[1].y(), TRV[1].z());
    septernion TR(septernion(CofG_ + TRV[0])*R*septernion(-CofG_));

    Info<< "solidBodyMotionFunctions::tabulated6DoFMotion::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::tabulated6DoFMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    // If the timeDataFileName has changed read the file

    fileName newTimeDataFileName
    (
        fileName(SBMFCoeffs_.lookup("timeDataFileName")).expand()
    );

    if (newTimeDataFileName != timeDataFileName_)
    {
        timeDataFileName_ = newTimeDataFileName;

        IFstream dataStream(timeDataFileName_);

        if (dataStream.good())
        {
            List<Tuple2<scalar, translationRotationVectors> > timeValues
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
            FatalErrorIn
            (
                "solidBodyMotionFunctions::tabulated6DoFMotion::read"
                "(const dictionary&)"
            )   << "Cannot open time data file " << timeDataFileName_
                << exit(FatalError);
        }
    }

    SBMFCoeffs_.lookup("CofG") >> CofG_;

    return true;
}


// ************************************************************************* //
