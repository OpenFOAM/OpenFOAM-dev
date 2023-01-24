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

#include "SDA.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(SDA, 0);
    addToRunTimeSelectionTable(solidBodyMotionFunction, SDA, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::SDA::SDA
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    CofG_(SBMFCoeffs_.lookup("CofG"))
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::SDA::~SDA()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion Foam::solidBodyMotionFunctions::SDA::transformation() const
{
    scalar time = time_.value();

    scalar Tpi = Tp_ + dTp_*(time/dTi_);   // Current roll period [sec]
    scalar wr = twoPi/Tpi; // Current Freq [/sec]

    // Current Phase for roll [rad]
    scalar r = dTp_/dTi_;
    scalar u = Tp_ + r*time;
    scalar phr = twoPi*((Tp_/u - 1) + log(mag(u)) - log(Tp_))/r;

    // Current Phase for Sway [rad]
    scalar phs = phr + pi;

    // Current Phase for Heave [rad]
    scalar phh = phr + piByTwo;

    scalar rollA = max(rollAmax_*exp(-sqr(Tpi - Tpn_)/(2*Q_)), rollAmin_);

    vector T
    (
        0,
        swayA_*(sin(wr*time + phs) - sin(phs)),
        heaveA_*(sin(wr*time + phh) - sin(phh))
    );
    quaternion R(quaternion::XYZ, vector(rollA*sin(wr*time + phr), 0, 0));
    septernion TR(septernion(-CofG_ - T)*R*septernion(CofG_));

    DebugInFunction << "Time = " << time << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::SDA::read(const dictionary& SBMFCoeffs)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("CofG") >> CofG_;
    SBMFCoeffs_.lookup("lambda") >> lambda_;
    SBMFCoeffs_.lookup("rollAmax") >> rollAmax_;
    SBMFCoeffs_.lookup("rollAmin") >> rollAmin_;
    SBMFCoeffs_.lookup("heaveA") >> heaveA_;
    SBMFCoeffs_.lookup("swayA") >> swayA_;
    SBMFCoeffs_.lookup("Q") >> Q_;
    SBMFCoeffs_.lookup("Tp") >> Tp_;
    SBMFCoeffs_.lookup("Tpn") >> Tpn_;
    SBMFCoeffs_.lookup("dTi") >> dTi_;
    SBMFCoeffs_.lookup("dTp") >> dTp_;

    // Rescale parameters according to the given scale parameter
    if (lambda_ > 1 + small)
    {
        heaveA_ /= lambda_;
        swayA_ /= lambda_;
        Tp_ /= sqrt(lambda_);
        Tpn_ /= sqrt(lambda_);
        dTi_ /= sqrt(lambda_);
        dTp_ /= sqrt(lambda_);
    }

    return true;
}


// ************************************************************************* //
