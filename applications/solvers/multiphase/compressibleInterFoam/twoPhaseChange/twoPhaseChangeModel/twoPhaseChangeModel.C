/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "twoPhaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseChangeModel, 0);
    defineRunTimeSelectionTable(twoPhaseChangeModel, dictionary);
}

const Foam::word Foam::twoPhaseChangeModel::phaseChangePropertiesName
(
    "phaseChangeProperties"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::twoPhaseChangeModel::createIOobject
(
    const compressibleTwoPhaseMixture& mixture
) const
{
    typeIOobject<IOdictionary> io
    (
        phaseChangePropertiesName,
        mixture.alpha1().mesh().time().constant(),
        mixture.alpha1().mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseChangeModel::twoPhaseChangeModel
(
    const word& type,
    const compressibleTwoPhaseMixture& mixture
)
:
    IOdictionary(createIOobject(mixture)),
    mixture_(mixture),
    twoPhaseChangeModelCoeffs_(optionalSubDict(type + "Coeffs"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::twoPhaseChangeModel::correct()
{}


bool Foam::twoPhaseChangeModel::read()
{
    if (regIOobject::read())
    {
        twoPhaseChangeModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
