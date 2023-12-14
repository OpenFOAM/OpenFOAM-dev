/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "splitCell.H"
#include "error.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from cell number and parent
Foam::splitCell::splitCell(const label celli, splitCell* parent)
:
    celli_(celli),
    parent_(parent),
    master_(nullptr),
    slave_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::splitCell::~splitCell()
{
    splitCell* myParent = parent();

    if (myParent)
    {
        // Make sure parent does not refer to me anymore.
        if (myParent->master() == this)
        {
            myParent->master() = nullptr;
        }
        else if (myParent->slave() == this)
        {
            myParent->slave() = nullptr;
        }
        else
        {
            FatalErrorInFunction
                << " parent's master or slave pointer" << endl
                << "Cell:" << cellLabel() << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::splitCell::isMaster() const
{
    splitCell* myParent = parent();

    if (!myParent)
    {
        FatalErrorInFunction
            << "Cell:" << cellLabel() << abort(FatalError);

        return false;
    }
    else if (myParent->master() == this)
    {
        return true;
    }
    else if (myParent->slave() == this)
    {
        return false;
    }
    else
    {
        FatalErrorInFunction
            << " parent's master or slave pointer" << endl
            << "Cell:" << cellLabel() << abort(FatalError);

        return false;
    }
}


bool Foam::splitCell::isUnrefined() const
{
    return !master() && !slave();
}


Foam::splitCell* Foam::splitCell::getOther() const
{
    splitCell* myParent = parent();

    if (!myParent)
    {
        FatalErrorInFunction
            << "Cell:" << cellLabel() << abort(FatalError);

        return nullptr;
    }
    else if (myParent->master() == this)
    {
        return myParent->slave();
    }
    else if (myParent->slave() == this)
    {
        return myParent->master();
    }
    else
    {
        FatalErrorInFunction
            << " parent's master or slave pointer" << endl
            << "Cell:" << cellLabel() << abort(FatalError);

        return nullptr;
    }
}


// ************************************************************************* //
