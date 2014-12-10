/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "regionModelFunctionObjectList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::regionModelFunctionObjectList::regionModelFunctionObjectList
(
    regionModel& owner
)
:
    PtrList<regionModelFunctionObject>(),
    owner_(owner),
    dict_(dictionary::null)
{}


Foam::regionModels::regionModelFunctionObjectList::regionModelFunctionObjectList
(
    regionModel& owner,
    const dictionary& dict,
    const bool readFields
)
:
    PtrList<regionModelFunctionObject>(),
    owner_(owner),
    dict_(dict)
{
    if (readFields)
    {
        wordList modelNames(dict.toc());

        Info<< "    Selecting region model functions" << endl;

        if (modelNames.size() > 0)
        {
            this->setSize(modelNames.size());

            forAll(modelNames, i)
            {
                const word& modelName = modelNames[i];

                this->set
                (
                    i,
                    regionModelFunctionObject::New
                    (
                        dict,
                        owner,
                        modelName
                    )
                );
            }
        }
        else
        {
            Info<< "    none" << endl;
        }
    }
}


Foam::regionModels::regionModelFunctionObjectList::regionModelFunctionObjectList
(
    const regionModelFunctionObjectList& cfol
)
:
    PtrList<regionModelFunctionObject>(cfol),
    owner_(cfol.owner_),
    dict_(cfol.dict_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::regionModels::regionModelFunctionObjectList::
~regionModelFunctionObjectList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionModels::regionModelFunctionObjectList::preEvolveRegion()
{
    forAll(*this, i)
    {
        this->operator[](i).preEvolveRegion();
    }
}


void Foam::regionModels::regionModelFunctionObjectList::postEvolveRegion()
{
    forAll(*this, i)
    {
        this->operator[](i).postEvolveRegion();
    }
}


// ************************************************************************* //
