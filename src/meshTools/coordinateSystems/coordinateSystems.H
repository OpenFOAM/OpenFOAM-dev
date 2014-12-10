/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::coordinateSystems

Description
    Provides a centralized coordinateSystem collection.

Note
    Mixing normal constructors and the coordinateSystems::New constructor
    may yield unexpected results.

    \verbatim
        1
        (
        cat1
        {
            coordinateSystem  system_10;
            porosity        0.781;
            Darcy
            {
                d   d [0 -2 0 0 0]  (-1000 -1000 0.50753e+08);
                f   f [0 -1 0 0 0]  (-1000 -1000 12.83);
            }
        }
        )
    \endverbatim

    For this to work correctly, the coordinateSystem constructor must be
    supplied with both a dictionary and an objectRegistry.

SourceFiles
    coordinateSystems.C

\*---------------------------------------------------------------------------*/
#ifndef coordinateSystems_H
#define coordinateSystems_H

#include "coordinateSystem.H"
#include "IOPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class coordinateSystems Declaration
\*---------------------------------------------------------------------------*/

class coordinateSystems
:
    public IOPtrList<coordinateSystem>
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        coordinateSystems(const coordinateSystems&);

        //- Disallow default bitwise assignment
        void operator=(const coordinateSystems&);


public:

    //- Runtime type information
    TypeName("coordinateSystems");

    // Constructors

        //- Read construct from IOobject
        explicit coordinateSystems(const IOobject&);

        //- Construct from IOobject and a PtrList
        coordinateSystems
        (
            const IOobject&,
            const PtrList<coordinateSystem>&
        );

        //- Construct from IOobject and transferring the PtrList contents
        coordinateSystems
        (
            const IOobject&,
            const Xfer<PtrList<coordinateSystem> >&
        );


    // Selectors

        //- Return previously registered or read construct from "constant"
        static const coordinateSystems& New(const objectRegistry&);


    // Member Functions

        //- Find and return indices for all matches
        labelList findIndices(const keyType& key) const;

        //- Find and return index for the first match, return -1 if not found
        label findIndex(const keyType& key) const;

        //- Find and return index for the first match, returns -1 if not found
        // \deprecated use findIndex() instead (deprecated Jul 2010)
        label find(const keyType& key) const;

        //- Search for given key
        bool found(const keyType& key) const;

        //- Return the table of contents (list of all keywords)
        wordList toc() const;

        //- write data
        bool writeData(Ostream&, bool subDict=true) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
