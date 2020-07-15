/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "diameterModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diameterModel, 0);
    defineRunTimeSelectionTable(diameterModel, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

Foam::volScalarField& Foam::diameterModel::dRef()
{
    if (!dPtr_.valid())
    {
        dPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("d", phase_.name()),
                    phase_.time().timeName(),
                    phase_.mesh()
                ),
                phase_.mesh(),
                dimensionedScalar(dimLength, 0)
            )
        );
    }

    return dPtr_();
}

Foam::volScalarField& Foam::diameterModel::aRef()
{
    if (!aPtr_.valid())
    {
        aPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("a", phase_.name()),
                    phase_.time().timeName(),
                    phase_.mesh()
                ),
                phase_.mesh(),
                dimensionedScalar(dimless/dimLength, 0)
            )
        );
    }

    return aPtr_();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModel::diameterModel
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    diameterProperties_(diameterProperties),
    phase_(phase),
    dPtr_(nullptr),
    aPtr_(nullptr)
{
    if (diameterProperties.lookupOrDefault("storeD", false))
    {
        dRef();
    }
    if (diameterProperties.lookupOrDefault("storeA", false))
    {
        aRef();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModel::~diameterModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::diameterModel::d() const
{
    if (dPtr_.valid())
    {
        return dPtr_();
    }
    else
    {
        return calcD();
    }
}


Foam::tmp<Foam::volScalarField> Foam::diameterModel::a() const
{
    if (aPtr_.valid())
    {
        return aPtr_();
    }
    else
    {
        return calcA();
    }
}


void Foam::diameterModel::correct()
{
    if (dPtr_.valid())
    {
        tmp<volScalarField> td = calcD();
        if (td.isTmp())
        {
            dPtr_() = td;
        }
    }
    if (aPtr_.valid())
    {
        tmp<volScalarField> tA = calcA();
        if (tA.isTmp())
        {
            aPtr_() = tA;
        }
    }
}


bool Foam::diameterModel::read(const dictionary& phaseProperties)
{
    diameterProperties_ = phaseProperties.optionalSubDict(type() + "Coeffs");

    return true;
}


// ************************************************************************* //
