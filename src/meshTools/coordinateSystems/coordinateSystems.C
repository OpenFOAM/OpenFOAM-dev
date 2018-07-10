/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "coordinateSystems.H"
#include "IOPtrList.H"
#include "Time.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coordinateSystems, 0);
    defineTemplateTypeNameAndDebug(IOPtrList<coordinateSystem>, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateSystems::coordinateSystems(const IOobject& io)
:
    IOPtrList<coordinateSystem>(io)
{}


Foam::coordinateSystems::coordinateSystems
(
    const IOobject& io,
    const PtrList<coordinateSystem>& lst
)
:
    IOPtrList<coordinateSystem>(io, lst)
{}


Foam::coordinateSystems::coordinateSystems
(
    const IOobject& io,
    const Xfer<PtrList<coordinateSystem>>& lst
)
:
    IOPtrList<coordinateSystem>(io, lst)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// Read construct from registry, or return previously registered
const Foam::coordinateSystems& Foam::coordinateSystems::New
(
    const objectRegistry& obr
)
{
    if (obr.foundObject<coordinateSystems>(typeName))
    {
        return obr.lookupObject<coordinateSystems>(typeName);
    }
    else
    {
        return obr.store
        (
            new coordinateSystems
            (
                IOobject
                (
                    typeName,
                    obr.time().constant(),
                    obr,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::coordinateSystems::findIndices(const keyType& key) const
{
    labelList indices;
    if (key.isPattern())
    {
        indices = findStrings(key, toc());
    }
    else
    {
        indices.setSize(size());
        label nFound = 0;
        forAll(*this, i)
        {
            if (key == operator[](i).name())
            {
                indices[nFound++] = i;
            }
        }
        indices.setSize(nFound);
    }

    return indices;
}


Foam::label Foam::coordinateSystems::findIndex(const keyType& key) const
{
    if (key.isPattern())
    {
        labelList indices = findIndices(key);
        // return first element
        if (!indices.empty())
        {
            return indices[0];
        }
    }
    else
    {
        forAll(*this, i)
        {
            if (key == operator[](i).name())
            {
                return i;
            }
        }
    }

    return -1;
}


bool Foam::coordinateSystems::found(const keyType& key) const
{
    return findIndex(key) != -1;
}


Foam::wordList Foam::coordinateSystems::toc() const
{
    wordList keywords(size());

    forAll(*this, i)
    {
        keywords[i] = operator[](i).name();
    }

    return keywords;
}


bool Foam::coordinateSystems::writeData(Ostream& os) const
{
    os << nl << size() << nl << token::BEGIN_LIST;

    forAll(*this, i)
    {
        os << nl;
        operator[](i).writeDict(os, true);
    }

    os << token::END_LIST << nl;

    return os.good();
}


// ************************************************************************* //
