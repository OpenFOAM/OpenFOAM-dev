/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "removeRegisteredObject.H"
#include "dictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(removeRegisteredObject, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::removeRegisteredObject::removeRegisteredObject
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    objectNames_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::removeRegisteredObject::~removeRegisteredObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::removeRegisteredObject::read(const dictionary& dict)
{
    dict.lookup("objectNames") >> objectNames_;
}


void Foam::removeRegisteredObject::execute()
{
    forAll(objectNames_, i)
    {
        if (obr_.foundObject<regIOobject>(objectNames_[i]))
        {
            const regIOobject& obj =
                obr_.lookupObject<regIOobject>(objectNames_[i]);

            if (obj.ownedByRegistry())
            {
                Info<< type() << " " << name_ << " output:" << nl
                    << "    removing object " << obj.name() << nl
                    << endl;

                const_cast<regIOobject&>(obj).release();
                delete &obj;
            }
        }
    }
}


void Foam::removeRegisteredObject::end()
{
    execute();
}


void Foam::removeRegisteredObject::timeSet()
{
    // Do nothing - only valid on execute
}


void Foam::removeRegisteredObject::write()
{
    // Do nothing - only valid on execute
}


// ************************************************************************* //
