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

Class
    Foam::CaCO3

Description
    Calcium carbonate (limestone)

SourceFiles
    CaCO3.C

\*---------------------------------------------------------------------------*/

#ifndef solid_CaCO3_H
#define solid_CaCO3_H

#include "solidProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class CaCO3;

Ostream& operator<<
(
    Ostream&,
    const CaCO3&
);


/*---------------------------------------------------------------------------*\
                          Class CaCO3 Declaration
\*---------------------------------------------------------------------------*/

class CaCO3
:
    public solidProperties
{

public:

    //- Runtime type information
    TypeName("CaCO3");


    // Constructors

        //- Construct null
        CaCO3();

        //- Construct from solidProperties
        CaCO3(const solidProperties& s);

        //- Construct from Istream
        CaCO3(Istream& is);

        //- Construct from dictionary
        CaCO3(const dictionary& dict);

        //- Construct copy
        CaCO3(const CaCO3& s);

        //- Construct and return clone
        virtual autoPtr<solidProperties> clone() const
        {
            return autoPtr<solidProperties>(new CaCO3(*this));
        }


    // I-O

        //- Write the function coefficients
        void writeData(Ostream& os) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const CaCO3& s);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
