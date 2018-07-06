/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "orderedPhasePair.H"
#include "aspectRatioModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::orderedPhasePair::orderedPhasePair
(
    const phaseModel& dispersed,
    const phaseModel& continuous,
    const dimensionedVector& g,
    const scalarTable& sigmaTable,
    const dictTable& aspectRatioTable
)
:
    phasePair(dispersed, continuous, g, sigmaTable, true)
{
    if (aspectRatioTable.found(*this))
    {
        aspectRatio_.set
        (
            aspectRatioModel::New
            (
                aspectRatioTable[*this],
                *this
            ).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::orderedPhasePair::~orderedPhasePair()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::phaseModel& Foam::orderedPhasePair::dispersed() const
{
    return phase1();
}


const Foam::phaseModel& Foam::orderedPhasePair::continuous() const
{
    return phase2();
}


Foam::word Foam::orderedPhasePair::name() const
{
    word namec(continuous().name());
    namec[0] = toupper(namec[0]);
    return dispersed().name() + "In" + namec;
}


Foam::tmp<Foam::volScalarField> Foam::orderedPhasePair::E() const
{
    if (!aspectRatio_.valid())
    {
        FatalErrorInFunction
            << "Aspect ratio model not specified for " << *this << "."
            << exit(FatalError);
    }

    return aspectRatio_->E();
}


// ************************************************************************* //
