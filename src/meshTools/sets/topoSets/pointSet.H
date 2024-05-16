/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::pointSet

Description
    A set of point labels.

SourceFiles
    pointSet.C

\*---------------------------------------------------------------------------*/

#ifndef pointSet_H
#define pointSet_H

#include "topoSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class pointSet Declaration
\*---------------------------------------------------------------------------*/

class pointSet
:
    public topoSet
{
public:

    //- Runtime type information
    TypeName("pointSet");


    // Constructors

        //- Construct from IOobject
        pointSet(const IOobject& obj);

        //- Construct from objectRegistry and name
        pointSet
        (
            const polyMesh& mesh,
            const word& name,
            readOption r=MUST_READ,
            writeOption w=NO_WRITE
        );

        //- Construct from additional size of labelHashSet
        pointSet
        (
            const polyMesh& mesh,
            const word& name,
            const label,
            writeOption w=NO_WRITE
        );

        //- Construct from additional labelHashSet
        pointSet
        (
            const polyMesh& mesh,
            const word& name,
            const topoSet&,
            writeOption w=NO_WRITE
        );

        //- Construct from additional labelHashSet
        pointSet
        (
            const polyMesh& mesh,
            const word& name,
            const labelHashSet&,
            writeOption w=NO_WRITE
        );


    //- Destructor
    virtual ~pointSet();


    // Member Functions

        //- Sync set across coupled patches. Adds coupled points to set.
        virtual void sync(const polyMesh& mesh);

        //- Return max index+1.
        virtual label maxSize(const polyMesh& mesh) const;

        //- Update any stored data for new labels
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Update any stored data for new labels
        // virtual void topoChange(const polyTopoChange& meshMod);

        //- Write maxLen items with label and coordinates.
        virtual void writeDebug
        (
            Ostream& os,
            const primitiveMesh&,
            const label maxLen
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
