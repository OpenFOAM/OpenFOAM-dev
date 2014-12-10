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

#include "fieldAverageItem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    const word fieldAverageItem::EXT_MEAN = "Mean";
    const word fieldAverageItem::EXT_PRIME2MEAN = "Prime2Mean";

    template<>
    const char* Foam::NamedEnum
    <
        Foam::fieldAverageItem::baseType,
        2
    >::names[] =
    {
        "iteration",
        "time"
    };
}


const Foam::NamedEnum<Foam::fieldAverageItem::baseType, 2>
    Foam::fieldAverageItem::baseTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldAverageItem::fieldAverageItem()
:
    active_(false),
    fieldName_("unknown"),
    mean_(0),
    meanFieldName_("unknown"),
    prime2Mean_(0),
    prime2MeanFieldName_("unknown"),
    base_(ITER),
    window_(-1.0),
    windowName_("")
{}


Foam::fieldAverageItem::fieldAverageItem(const fieldAverageItem& faItem)
:
    active_(faItem.active_),
    fieldName_(faItem.fieldName_),
    mean_(faItem.mean_),
    meanFieldName_(faItem.meanFieldName_),
    prime2Mean_(faItem.prime2Mean_),
    prime2MeanFieldName_(faItem.prime2MeanFieldName_),
    base_(faItem.base_),
    window_(faItem.window_),
    windowName_(faItem.windowName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldAverageItem::~fieldAverageItem()
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::fieldAverageItem::operator=(const fieldAverageItem& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "Foam::fieldAverageItem::operator=(const Foam::fieldAverageItem&)"
        )   << "Attempted assignment to self" << nl
            << abort(FatalError);
    }

    // Set updated values
    active_ = rhs.active_;
    fieldName_ = rhs.fieldName_;
    mean_ = rhs.mean_;
    meanFieldName_ = rhs.meanFieldName_;
    prime2Mean_ = rhs.prime2Mean_;
    prime2MeanFieldName_ = rhs.prime2MeanFieldName_;
    base_ = rhs.base_;
    window_ = rhs.window_;
    windowName_ = rhs.windowName_;
}


// ************************************************************************* //
