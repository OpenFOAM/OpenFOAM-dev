/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"

#include "scalar.H"

#include "IOstreams.H"
#include "Dictionary.H"
#include "PtrDictionary.H"

using namespace Foam;

class ent
:
    public Dictionary<ent>::link
{
    word keyword_;
    int i_;

public:

    ent(const word& keyword, int i)
    :
        keyword_(keyword),
        i_(i)
    {}

    const word& keyword() const
    {
        return keyword_;
    }

    friend Ostream& operator<<(Ostream& os, const ent& e)
    {
        os  << e.keyword_ << ' ' << e.i_ << endl;
        return os;
    }
};


class Scalar
{
    scalar data_;

public:

    Scalar()
    :
        data_(0)
    {}

    Scalar(scalar val)
    :
        data_(val)
    {}

    ~Scalar()
    {
        Info<<"delete Scalar: " << data_ << endl;
    }

    friend Ostream& operator<<(Ostream& os, const Scalar& val)
    {
        os  << val.data_;
        return os;
    }
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    Dictionary<ent>* dictPtr = new Dictionary<ent>;
    Dictionary<ent>& dict = *dictPtr;

    for (int i = 0; i<10; i++)
    {
        ent* ePtr = new ent(word("ent") + name(i), i);
        dict.append(ePtr->keyword(), ePtr);
        dict.swapUp(ePtr);
    }

    Info<< dict << endl;

    dict.swapDown(dict.first());

    forAllConstIter(Dictionary<ent>, dict, iter)
    {
        Info<< "element : " << *iter;
    }

    Info<< "keys: " << dict.toc() << endl;

    delete dictPtr;

    Dictionary<ent> dict2;

    for (int i = 0; i<10; i++)
    {
        ent* ePtr = new ent(word("ent") + name(i), i);
        dict2.append(ePtr->keyword(), ePtr);
        dict2.swapUp(ePtr);
    }

    Info<< "dict:\n" << dict2 << endl;

    Info<< nl << "Testing transfer: " << nl << endl;
    Info<< "original: " << dict2 << endl;

    Dictionary<ent> newDict;
    newDict.transfer(dict2);

    Info<< nl << "source: " << dict2 << nl
        << "keys: " << dict2.toc() << nl
        << "target: " << newDict << nl
        << "keys: " << newDict.toc() << endl;


    PtrDictionary<Scalar> scalarDict;
    for (int i = 0; i<10; i++)
    {
        word key("ent" + name(i));
        scalarDict.insert(key, new Scalar(1.3*i));
    }

    Info<< nl << "scalarDict1: " << endl;
    forAllConstIter(PtrDictionary<Scalar>, scalarDict, iter)
    {
        Info<< " = " << iter() << endl;
    }

    PtrDictionary<Scalar> scalarDict2;
    for (int i = 8; i<15; i++)
    {
        word key("ent" + name(i));
        scalarDict2.insert(key, new Scalar(1.3*i));
    }
    Info<< nl << "scalarDict2: " << endl;
    forAllConstIter(PtrDictionary<Scalar>, scalarDict2, iter)
    {
        Info<< "elem = " << *iter << endl;
    }

    scalarDict.transfer(scalarDict2);


    Scalar* p = scalarDict.lookupPtr("ent8");

    // This does not (yet) work
    // Scalar* q = scalarDict.remove("ent10");

    if (p)
    {
        Info<< "found: " << *p << endl;
    }
    else
    {
        Info<< "no p: " << endl;
    }

    scalarDict.clear();

    // Info<< " = " << *iter << endl;



    Info<< nl << "Done." << endl;
    return 0;
}


// ************************************************************************* //
