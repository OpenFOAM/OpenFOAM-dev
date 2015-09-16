/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
#include "PtrListDictionary.H"

using namespace Foam;

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
    PtrListDictionary<Scalar> scalarDict(10);
    forAll(scalarDict, i)
    {
        word key("ent" + name(i));
        scalarDict.set(i, key, new Scalar(1.3*i));
    }

    Info<< nl << "scalarDict1: " << endl;
    forAll(scalarDict, i)
    {
        Info<< "elem " << i << " = " << scalarDict[i] << endl;
    }

    Scalar* ent8Ptr = scalarDict.lookupPtr("ent8");

    Info<< "ent8 = " << *ent8Ptr << endl;

    PtrListDictionary<Scalar> scalarDict2(15);
    forAll(scalarDict2, i)
    {
        word key("ent" + name(i));
        scalarDict2.set(i, key, new Scalar(1.3*i));
    }
    Info<< nl << "scalarDict2: " << endl;
    forAll(scalarDict2, i)
    {
        Info<< "elem " << i << " = " << scalarDict2[i] << endl;
    }

    scalarDict.transfer(scalarDict2);

    Scalar* p = scalarDict.lookupPtr("ent8");

    if (p)
    {
        Info<< "found: " << *p << endl;
    }
    else
    {
        Info<< "no p: " << endl;
    }

    scalarDict.clear();

    Info<< nl << "Done." << endl;
    return 0;
}


// ************************************************************************* //
