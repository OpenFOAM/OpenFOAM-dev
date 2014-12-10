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

#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionary* dimensionedConstantsPtr_(NULL);

dictionary& dimensionedConstants()
{
    return debug::switchSet
    (
        "DimensionedConstants",
        dimensionedConstantsPtr_
    );
}


dimensionedScalar dimensionedConstant
(
    const word& group,
    const word& varName
)
{
    dictionary& dict = dimensionedConstants();

    // Check that the entries exist.
    // Note: should make FatalError robust instead!

    if (!dict.found("unitSet"))
    {
        std::cerr<< "Cannot find unitSet in dictionary " << dict.name()
            << std::endl;
    }

    const word unitSetCoeffs(word(dict.lookup("unitSet")) + "Coeffs");

    if (!dict.found(unitSetCoeffs))
    {
        std::cerr<< "Cannot find " << unitSetCoeffs << " in dictionary "
            << dict.name() << std::endl;
    }

    dictionary& unitDict = dict.subDict(unitSetCoeffs);

    dictionary& groupDict = unitDict.subDict(group);

    return dimensionedScalar(groupDict.lookup(varName));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
