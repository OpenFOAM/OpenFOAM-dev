/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "parcelThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parcelThermo::parcelThermo(const fluidThermo& carrierThermo)
:
    liquids_(nullptr),
    solids_(nullptr)
{
    Info<< "Creating component thermo properties:" << endl;

    if (carrierThermo.properties().found("liquids"))
    {
        liquids_ = liquidMixtureProperties::New
        (
            carrierThermo.properties().subDict("liquids")
        );
        Info<< "    liquids - " << liquids_->components().size()
            << " components" << endl;
    }
    else
    {
        Info<< "    no liquid components" << endl;
    }

    if (carrierThermo.properties().found("solids"))
    {
        solids_ = solidMixtureProperties::New
        (
            carrierThermo.properties().subDict("solids")
        );
        Info<< "    solids - " << solids_->components().size()
            << " components" << endl;
    }
    else
    {
        Info<< "    no solid components" << endl;
    }
}


Foam::parcelThermo::parcelThermo(const parcelThermo& pThermo)
:
    liquids_(pThermo.liquids_, false),
    solids_(pThermo.solids_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::parcelThermo::~parcelThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::liquidMixtureProperties& Foam::parcelThermo::liquids() const
{
    if (!liquids_.valid())
    {
        FatalErrorInFunction
            << "liquids requested, but object is not allocated"
            << abort(FatalError);
    }

    return liquids_();
}


const Foam::solidMixtureProperties& Foam::parcelThermo::solids() const
{
    if (!solids_.valid())
    {
        FatalErrorInFunction
            << "solids requested, but object is not allocated"
            << abort(FatalError);
    }

    return solids_();
}


Foam::label Foam::parcelThermo::liquidId
(
    const word& cmptName,
    bool allowNotfound
) const
{
    forAll(liquids().components(), i)
    {
        if (cmptName == liquids_->components()[i])
        {
            return i;
        }
    }

    if (!allowNotfound)
    {
        FatalErrorInFunction
            << "Unknown liquid component " << cmptName << ". Valid liquids are:"
            << nl << liquids_->components() << exit(FatalError);
    }

    return -1;
}


Foam::label Foam::parcelThermo::solidId
(
    const word& cmptName,
    bool allowNotfound
) const
{
    forAll(solids().components(), i)
    {
        if (cmptName == solids_->components()[i])
        {
            return i;
        }
    }

    if (!allowNotfound)
    {
        FatalErrorInFunction
            << "Unknown solid component " << cmptName << ". Valid solids are:"
            << nl << solids_->components() << exit(FatalError);
    }

    return -1;
}


// ************************************************************************* //
