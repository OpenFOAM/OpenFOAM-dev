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
    Foam::zoneGeneratorList

Description
    List of zoneGenerators

    Used by Foam::zonesGenerator to generate and store the zones and also by the
    \c union, \c difference and \c intersection \c zoneGenerators to construct
    the list of zones they combine.

See also
    Foam::zonesGenerator
    Foam::zoneGenerator
    Foam::zoneSet

SourceFiles
    zoneGeneratorList.C

\*---------------------------------------------------------------------------*/
#ifndef zoneGeneratorList_H
#define zoneGeneratorList_H

#include "zoneGenerator.H"
#include "PtrListDictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class zoneGeneratorList Declaration
\*---------------------------------------------------------------------------*/

class zoneGeneratorList
:
    public PtrListDictionary<zoneGenerator>
{
protected:

    // Protected Data

        //- Reference to the polyMesh
        const polyMesh& mesh_;


    // Protected Member Functions

        void read(const dictionary& dict);


public:

    // Constructors

        //- Construct from mesh
        zoneGeneratorList(const polyMesh& mesh);

        //- Construct from mesh and dictionary
        zoneGeneratorList(const polyMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        zoneGeneratorList(const zoneGeneratorList&) = delete;


    // Member Functions

        //- Return true if any of the zoneGenerators updates any of the zones
        // following mesh point motion
        bool moveUpdate() const;

        //- Generate and register the zones
        void generate();

        //- Update the zones following mesh point movement
        void movePoints();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const zoneGeneratorList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
