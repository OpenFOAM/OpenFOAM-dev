/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::surfacePatch

Description
    'Patch' on surface as subset of triSurface.

SourceFiles
    surfacePatch.C

\*---------------------------------------------------------------------------*/

#ifndef surfacePatch_H
#define surfacePatch_H

#include "geometricSurfacePatch.H"
#include "word.H"
#include "label.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class surfacePatch;

Ostream& operator<<(Ostream&, const surfacePatch&);


/*---------------------------------------------------------------------------*\
                           Class surfacePatch Declaration
\*---------------------------------------------------------------------------*/

class surfacePatch
:
    public geometricSurfacePatch
{
    // Private Data

        //- Size of this patch in the triSurface face list
        label size_;

        //- Start label of this patch in the triSurface face list
        label start_;


public:

    //- Runtime type information
    ClassName("surfacePatch");


    // Constructors

        //- Construct null
        surfacePatch();

        //- Construct from components
        surfacePatch
        (
            const word& geometricType,
            const word& name,
            const label size,
            const label start,
            const label index
        );

        //- Construct from Istream
        surfacePatch(Istream& is, const label index);

        //- Construct from dictionary
        surfacePatch
        (
            const word& name,
            const dictionary& dict,
            const label index
        );


    // Member Functions

        //- Return start label of this patch in the polyMesh face list
        label start() const
        {
            return start_;
        }

        //- Return start label of this patch in the polyMesh face list
        label& start()
        {
            return start_;
        }

        //- Return size of this patch in the polyMesh face list
        label size() const
        {
            return size_;
        }

        //- Return size of this patch in the polyMesh face list
        label& size()
        {
            return size_;
        }

        //- Write
        void write(Ostream&) const;

        //- Write dictionary
        void writeDict(Ostream&) const;


    // Member Operators

        bool operator!=(const surfacePatch&) const;

        //- compare.
        bool operator==(const surfacePatch&) const;

    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const surfacePatch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
