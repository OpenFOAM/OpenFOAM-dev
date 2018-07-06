/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::functionObjects::volRegion

Description
    Volume (cell) region selection class.

    Examples of function object specification:
    \verbatim
    volRegion0
    {
        .
        .
        regionType      cellZone;
        name            c0;
        .
        .
    }

    volRegionAll
    {
        .
        .
        regionType      all;
        .
        .
    }
    \endverbatim

Usage
    \table
        Property     | Description                | Required     | Default value
        regionType   | cellZone or all              | no | all
        name         | Name of cellZone if required | no |
    \endtable

See also
    Foam::functionObject

SourceFiles
    volRegion.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_volRegion_H
#define functionObjects_volRegion_H

#include "writeFile.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class volRegion Declaration
\*---------------------------------------------------------------------------*/

class volRegion
{
    // Private member data

        const fvMesh& mesh_;

        // Cache integral properties of the region for writeFileHeader
        label nCells_;
        scalar V_;


public:

    // Public data types

        //- Region type enumeration
        enum regionTypes
        {
            vrtCellZone,
            vrtAll
        };

        //- Region type names
        static const NamedEnum<regionTypes, 2> regionTypeNames_;


protected:

    // Protected data

        //- Region type
        regionTypes regionType_;

        //- Region name (patch, zone, etc.)
        word regionName_;

        //- Region ID (patch ID, zone ID, etc.)
        label regionID_;


    // Protected Member Functions

        //- Output file header information
        void writeFileHeader(const writeFile& wf, Ostream& file);


public:

    //- Run-time type information
    TypeName("volRegion");


    // Constructors

        //- Construct from fvMesh and dictionary
        volRegion
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~volRegion();


    // Public Member Functions

        //- Read from dictionary
        bool read(const dictionary&);

        //- Return the region type
        inline const regionTypes& regionType() const;

        //- Return the local list of cell IDs
        const labelList& cellIDs() const;

        //- Return the number of cells in the region
        label nCells() const;

        //- Return total volume of the region
        scalar V() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volRegionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
