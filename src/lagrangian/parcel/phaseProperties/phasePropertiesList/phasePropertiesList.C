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

#include "phasePropertiesList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phasePropertiesList::phasePropertiesList()
:
    props_(),
    phaseTypeNames_(),
    stateLabels_()
{}


Foam::phasePropertiesList::phasePropertiesList
(
    Istream& is,
    const wordList& gasNames,
    const wordList& liquidNames,
    const wordList& solidNames
)
:
    props_(is),
    phaseTypeNames_(),
    stateLabels_()
{
    forAll(props_, i)
    {
        props_[i].reorder(gasNames, liquidNames, solidNames);
    }

    phaseTypeNames_.setSize(props_.size());
    stateLabels_.setSize(props_.size());
    forAll(props_, i)
    {
        phaseTypeNames_[i] = props_[i].phaseTypeName();
        stateLabels_[i] = props_[i].stateLabel();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phasePropertiesList::~phasePropertiesList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::List<Foam::phaseProperties>&
Foam::phasePropertiesList::props() const
{
    return props_;
}


const Foam::wordList& Foam::phasePropertiesList::phaseTypes() const
{
    return phaseTypeNames_;
}


const Foam::wordList& Foam::phasePropertiesList::stateLabels() const
{
    return stateLabels_;
}


Foam::label Foam::phasePropertiesList::size() const
{
    return props_.size();
}


const Foam::phaseProperties&
Foam::phasePropertiesList::operator[](const label phaseI) const
{
    return props_[phaseI];
}


// ************************************************************************* //
