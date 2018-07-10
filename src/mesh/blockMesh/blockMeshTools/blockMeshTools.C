/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "blockMeshTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockMeshTools::read
(
    Istream& is,
    label& val,
    const dictionary& dict
)
{
    token t(is);
    if (t.isLabel())
    {
        val = t.labelToken();
    }
    else if (t.isWord())
    {
        const word& varName = t.wordToken();
        const entry* ePtr = dict.lookupScopedEntryPtr
        (
            varName,
            true,
            true
        );
        if (ePtr)
        {
            // Read as label
            val = Foam::readLabel(ePtr->stream());
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "Undefined variable "
                << varName << ". Valid variables are " << dict
                << exit(FatalIOError);
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Illegal token " << t.info()
            << " when trying to read label"
            << exit(FatalIOError);
    }

    is.fatalCheck
    (
        "operator>>(Istream&, List<T>&) : reading entry"
    );
}


Foam::label Foam::blockMeshTools::read
(
    Istream& is,
    const dictionary& dict
)
{
    label val;
    read(is, val, dict);
    return val;
}


void Foam::blockMeshTools::write
(
    Ostream& os,
    const label val,
    const dictionary& dict
)
{
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isStream())
        {
            label keyVal(Foam::readLabel(iter().stream()));
            if (keyVal == val)
            {
                os << iter().keyword();
                return;
            }
        }
    }
    os << val;
}


const Foam::keyType& Foam::blockMeshTools::findEntry
(
    const dictionary& dict,
    const label val
)
{
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isStream())
        {
            label keyVal(Foam::readLabel(iter().stream()));
            if (keyVal == val)
            {
                return iter().keyword();
            }
        }
    }

    return keyType::null;
}


// ************************************************************************* //
