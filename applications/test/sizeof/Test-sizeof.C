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

Description
    Test the sizeof various classes.

\*---------------------------------------------------------------------------*/

#include "bool.H"
#include "Switch.H"
#include "string.H"
#include "dictionary.H"
#include "nil.H"
#include "IOstreams.H"
#include "IStringStream.H"

namespace Foam
{
   class hasBoolClass
   {
   public:
      bool b_;

      hasBoolClass(const bool val=false)
      :
          b_(false)
      {}
   };

}


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    cout<<"sizeof\n------\n";
    {
        nil x;
        cout<<"nil:" << sizeof(x) << nl;
    }
    {
        bool x(0);
        cout<<"bool:" << sizeof(x) << nl;
    }
    {
        hasBoolClass x(true);
        cout<<"hasBoolClass:" << sizeof(x) << nl;
    }

    {
        Switch x("n");
        cout<<"Switch:" << sizeof(x) << nl;
        cout<<"Switch::switchType=" << sizeof(Switch::switchType) << nl;
    }

    {
        scalar x(0);
        cout<<"scalar:" << sizeof(x) << nl;
    }

    {
        label x(0);
        cout<<"label:" << sizeof(x) << nl;
    }

    {
        cout<<"int:" << sizeof(int) << nl;
        cout<<"long:" << sizeof(long) << nl;
        cout<<"float:" << sizeof(float) << nl;
        cout<<"double:" << sizeof(double) << nl;
    }


    Info << "---\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
