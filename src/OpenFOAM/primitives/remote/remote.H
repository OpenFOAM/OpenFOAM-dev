/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::remote

Description
    Struct for keeping processor, element (cell, face, point) index.

SourceFiles
    remoteI.H

\*---------------------------------------------------------------------------*/

#ifndef remote_H
#define remote_H

#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class remote;
inline bool operator==(const remote&, const remote&);
inline bool operator!=(const remote&, const remote&);
inline Istream& operator>>(Istream&, remote&);
inline Ostream& operator<<(Ostream&, const remote&);

/*---------------------------------------------------------------------------*\
                          Class remote Declaration
\*---------------------------------------------------------------------------*/

class remote
{
public:

    // Public Data

        //- Processor index
        label proci;

        //- Element index
        label elementi;


    // Constructors

        //- Construct null
        inline remote();

        //- Construct from components
        inline remote(const label, const label);

        //- Construct from stream
        inline remote(Istream& is);


    // Friend Operators

        //- Equality comparison
        friend bool operator==(const remote& a, const remote& b);

        //- Inequality comparison
        friend bool operator!=(const remote& a, const remote& b);


    // IOstream Operators

        //- Write to stream
        friend Ostream& operator<<(Ostream& os, const remote& p);

        //- Read from stream
        friend Istream& operator>>(Istream& is, remote& p);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
class pTraits<remote>
{
public:
    typedef remote cmptType;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "remoteI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
