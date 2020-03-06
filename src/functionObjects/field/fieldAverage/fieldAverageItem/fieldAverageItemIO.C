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

#include "fieldAverageItem.H"
#include "dictionaryEntry.H"
#include "IOobject.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldAverageItem::fieldAverageItem(Istream& is)
:
    fieldName_("unknown"),
    mean_(0),
    meanFieldName_("unknown"),
    prime2Mean_(0),
    prime2MeanFieldName_("unknown"),
    base_(baseType::iter),
    window_(-1.0)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::functionObjects::operator>>
(
    Istream& is,
    fieldAverageItem& faItem
)
{
    is.check
    (
        "Foam::Istream& Foam::operator>>"
        "(Foam::Istream&, Foam::functionObjects::fieldAverageItem&)"
    );

    token fieldNameToken(is);
    faItem.fieldName_ = fieldNameToken.wordToken();

    token nextToken(is);
    is.putBack(nextToken);

    if (nextToken.isPunctuation() && nextToken.pToken() == token::BEGIN_BLOCK)
    {
        const dictionary entry(dictionary::null, is);

        faItem.mean_ = entry.lookupOrDefault<Switch>("mean", true);
        faItem.prime2Mean_ = entry.lookupOrDefault<Switch>("prime2Mean", false);
        faItem.base_ =
            faItem.baseTypeNames_[entry.lookupOrDefault<word>("base", "time")];
        faItem.window_ = entry.lookupOrDefault<scalar>("window", -1);
        faItem.windowName_ = entry.lookupOrDefault<word>("windowName", "");
    }
    else
    {
        faItem.mean_ = true;
        faItem.prime2Mean_ = false;
        faItem.base_ = faItem.baseTypeNames_["time"];
        faItem.window_ = -1;
        faItem.windowName_ = "";
    }

    faItem.meanFieldName_ = IOobject::groupName
    (
        IOobject::member(faItem.fieldName_) + fieldAverageItem::meanExt,
        IOobject::group(faItem.fieldName_)
    );

    faItem.prime2MeanFieldName_ = IOobject::groupName
    (
        IOobject::member(faItem.fieldName_) + fieldAverageItem::prime2MeanExt,
        IOobject::group(faItem.fieldName_)
    );

    if ((faItem.window_ > 0) && (faItem.windowName_ != ""))
    {
        faItem.meanFieldName_ =
            faItem.meanFieldName_ + "_" + faItem.windowName_;

        faItem.prime2MeanFieldName_ =
            faItem.prime2MeanFieldName_ + "_" + faItem.windowName_;
    }

    return is;
}


// ************************************************************************* //
