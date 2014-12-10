/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "fieldValue.H"
#include "fvMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fieldValue, 0);
    defineRunTimeSelectionTable(fieldValue, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fieldValue::read(const dictionary& dict)
{
    if (active_)
    {
        dict_ = dict;

        log_ = dict.lookupOrDefault<Switch>("log", true);
        dict.lookup("fields") >> fields_;
        dict.lookup("valueOutput") >> valueOutput_;
    }
}


void Foam::fieldValue::write()
{
    if (active_)
    {
        functionObjectFile::write();

        Info(log_)<< type() << " " << name_ << " output:" << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldValue::fieldValue
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const word& valueType,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name, valueType),
    name_(name),
    obr_(obr),
    dict_(dict),
    active_(true),
    log_(true),
    sourceName_(dict.lookupOrDefault<word>("sourceName", "sampledSurface")),
    fields_(dict.lookup("fields")),
    valueOutput_(dict.lookup("valueOutput")),
    resultDict_(fileName("name"), dictionary::null)
{
    // Only active if obr is an fvMesh
    if (isA<fvMesh>(obr_))
    {
        read(dict);
    }
    else
    {
        WarningIn
        (
            "fieldValue::fieldValue"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name << nl
            << endl;
        active_ = false;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldValue::~fieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldValue::execute()
{
    // Do nothing
}


void Foam::fieldValue::end()
{
    // Do nothing
}


void Foam::fieldValue::timeSet()
{
    // Do nothing
}


void Foam::fieldValue::updateMesh(const mapPolyMesh&)
{
    // Do nothing
}


void Foam::fieldValue::movePoints(const polyMesh&)
{
    // Do nothing
}


// ************************************************************************* //
