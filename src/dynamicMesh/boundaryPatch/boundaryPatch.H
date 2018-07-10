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

Class
    Foam::boundaryPatch

Description
    Like polyPatch but without reference to mesh. patchIdentifier::index
    is not used. Used in boundaryMesh to hold data on patches.

SourceFiles
    boundaryPatch.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryPatch_H
#define boundaryPatch_H

#include "patchIdentifier.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class boundaryPatch;

Ostream& operator<<(Ostream&, const boundaryPatch&);


/*---------------------------------------------------------------------------*\
                           Class boundaryPatch Declaration
\*---------------------------------------------------------------------------*/

class boundaryPatch
:
    public patchIdentifier
{
    // Private data

        label size_;
        label start_;

public:

    // Constructors

        //- Construct from components
        boundaryPatch
        (
            const word& name,
            const label index,
            const label size,
            const label start,
            const word& physicalType = word::null
        );

        //- Construct from dictionary
        boundaryPatch
        (
            const word& name,
            const dictionary& dict,
            const label index
        );

        //- Construct as copy
        boundaryPatch(const boundaryPatch&);

        //- Construct as copy, resetting the index
        boundaryPatch(const boundaryPatch&, const label index);

        //- Clone
        autoPtr<boundaryPatch> clone() const;


    //- Destructor
    ~boundaryPatch();


    // Member Functions

        label size() const
        {
            return size_;
        }

        label& size()
        {
            return size_;
        }

        label start() const
        {
            return start_;
        }

        label& start()
        {
            return start_;
        }


        //- Write dictionary
        virtual void write(Ostream&) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const boundaryPatch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
