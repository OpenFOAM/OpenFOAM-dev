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
    Foam::zoneGenerator

Description
    Abstract base class for all zoneGenerators, providing runtime selection

See also
    Foam::zoneSet
    Foam::zoneGeneratorList

SourceFiles
    zoneGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef zoneGenerator_H
#define zoneGenerator_H

#include "zoneSet.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;

/*---------------------------------------------------------------------------*\
                          Class zoneGenerator Declaration
\*---------------------------------------------------------------------------*/

class zoneGenerator
{
protected:

    // Protected Data

        //- Name of zone generator
        const word name_;

        //- Dictionary cached for error context
        const dictionary dict_;

        //- Name of the zone (defaults to the name of the generator)
        const word zoneName_;

        //- Reference to the polyMesh
        const polyMesh& mesh_;

        //- Switch to update the zones if the mesh points are moved
        mutable bool moveUpdate_;


    // Protected Member Functions

        template<class ZoneType>
        labelList select
        (
            const ZoneType& zone,
            const vectorField& pts
        ) const;


public:

    friend class zoneGeneratorList;


    // Public Data Types

        //- Enumeration defining the file checking options
        enum zoneTypes
        {
            pointZoneType,
            cellZoneType,
            faceZoneType
        };

        static const NamedEnum<zoneTypes, 3> zoneTypesNames;


    //- Runtime type information
    TypeName("zoneGenerator");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            zoneGenerator,
            dictionary,
            (
                const word& name,
                const polyMesh& mesh,
                const dictionary& dict
            ),
            (name, mesh, dict)
        );


    // Constructors

        //- Construct from name, polyMesh and dictionary
        zoneGenerator
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        zoneGenerator(const zoneGenerator&) = delete;


    // Selectors

        //- Select constructed from name, mesh and dictionary
        static autoPtr<zoneGenerator> New
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Select constructed from mesh and first valid dictionary
        static autoPtr<zoneGenerator> New
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~zoneGenerator();


    // Member Functions

        //- Return the list of selected indices
        static labelList indices(const boolList& selected);

        //- Return true if the zoneGenerator updates any of the zones
        // following mesh point motion
        bool moveUpdate() const
        {
            return moveUpdate_;
        }

        //- Generate and return the zoneSet
        virtual zoneSet generate() const = 0;

        //- Regenerate the zoneSet following mesh point motion
        virtual zoneSet movePoints() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const zoneGenerator&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "zoneGeneratorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
