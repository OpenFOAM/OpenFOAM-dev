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

Application
    callBackTest

Description

\*---------------------------------------------------------------------------*/

#include "Callback.H"

using namespace Foam;

class callback
:
    public Callback<callback>
{
public:

    callback(CallbackRegistry<callback>& cbr)
    :
        Callback<callback>(cbr)
    {}

    ~callback()
    {}

    virtual const word& name() const
    {
        return word::null;
    }

    void testCallbackFunction() const
    {
        Info<< "calling testCallbackFunction for object " << name() << endl;
    }
};


class callbackRegistry
:
    public CallbackRegistry<callback>
{
public:

    callbackRegistry()
    {}

    ~callbackRegistry()
    {}

    void testCallbackFunction() const
    {
        forAllConstIter(callbackRegistry, *this, iter)
        {
            iter().testCallbackFunction();
        }
    }
};


class objectWithCallback
:
    public callback
{
    word name_;

public:

    objectWithCallback(const word& n, callbackRegistry& cbr)
    :
        callback(cbr),
        name_(n)
    {}

    virtual const word& name() const
    {
        return name_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    callbackRegistry cbr;

    objectWithCallback ob1("ob1", cbr);
    objectWithCallback ob2("ob2", cbr);

    cbr.testCallbackFunction();

    {
        objectWithCallback ob1("ob1", cbr);
        cbr.testCallbackFunction();
    }

    cbr.testCallbackFunction();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
