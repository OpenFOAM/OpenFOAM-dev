/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::generatedZoneSet

Description


See also
    Foam::zoneSet
    Foam::zoneGenerator

SourceFiles
    generatedZoneSet.C

\*---------------------------------------------------------------------------*/

#ifndef generatedZoneSet_H
#define generatedZoneSet_H

#include "zoneGenerator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class generatedZoneSet Declaration
\*---------------------------------------------------------------------------*/

class generatedZoneSet
:
    public zoneSet
{
    // Private Member Data

        autoPtr<zoneGenerator> zoneGenerator_;


public:

    // Constructors

        //- Construct null
        generatedZoneSet()
        {}

        //- Construct from name, mesh and dictionary
        generatedZoneSet
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from name, zoneType, mesh and dictionary
        generatedZoneSet
        (
            const word& name,
            const zoneGenerator::zoneTypes& zoneType,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        generatedZoneSet(const generatedZoneSet&) = delete;


    // Member Functions

        //- Read from dictionary
        void read
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Read from dictionary
        void read
        (
            const word& name,
            const zoneGenerator::zoneTypes& zoneType,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Set the zoneGenerator and generate the zoneSet
        void set(const autoPtr<zoneGenerator>& zg);

        //- Update for mesh motion
        bool movePoints();

        //- Redistribute or update using the given distribution map
        void distribute(const polyDistributionMap& map);

        //- Update topology using the given map
        void topoChange(const polyTopoChangeMap& map);

        //- Update from another mesh using the given map
        void mapMesh(const polyMeshMap& map);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const generatedZoneSet&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
