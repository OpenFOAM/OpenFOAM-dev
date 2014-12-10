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

#include "seriesProfile.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(seriesProfile, 0);
    addToRunTimeSelectionTable(profileModel, seriesProfile, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::seriesProfile::evaluateDrag
(
    const scalar& xIn,
    const List<scalar>& values
) const
{
    scalar result = 0.0;

    forAll(values, i)
    {
        result += values[i]*cos(i*xIn);
    }

    return result;
}


Foam::scalar Foam::seriesProfile::evaluateLift
(
    const scalar& xIn,
    const List<scalar>& values
) const
{
    scalar result = 0.0;

    forAll(values, i)
    {
        // note: first contribution always zero since sin(0) = 0, but
        // keep zero base to be consitent with drag coeffs
        result += values[i]*sin(i*xIn);
    }

    return result;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::seriesProfile::seriesProfile
(
    const dictionary& dict,
    const word& modelName
)
:
    profileModel(dict, modelName),
    CdCoeffs_(),
    ClCoeffs_()
{
    if (readFromFile())
    {
        IFstream is(fName_);
        is  >> CdCoeffs_ >> ClCoeffs_;
    }
    else
    {
        dict.lookup("CdCoeffs") >> CdCoeffs_;
        dict.lookup("ClCoeffs") >> ClCoeffs_;
    }


    if (CdCoeffs_.empty())
    {
        FatalErrorIn
        (
            "Foam::seriesProfile::seriesProfile"
            "("
                "const dictionary&, "
                "const word&"
            ")"
        )   << "CdCoeffs must be specified" << exit(FatalError);
    }
    if (ClCoeffs_.empty())
    {
        FatalErrorIn
        (
            "Foam::seriesProfile::seriesProfile"
            "("
                "const dictionary&, "
                "const word&"
            ")"
        )   << "ClCoeffs must be specified" << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::seriesProfile::Cdl(const scalar alpha, scalar& Cd, scalar& Cl) const
{
    Cd = evaluateDrag(alpha, CdCoeffs_);
    Cl = evaluateLift(alpha, ClCoeffs_);
}


// ************************************************************************* //
