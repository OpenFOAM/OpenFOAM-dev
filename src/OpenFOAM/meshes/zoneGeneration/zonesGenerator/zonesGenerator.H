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
    Foam::zonesGenerator

Description
    MeshObject zoneGeneratorList instantiated by polyMesh

    to provide optional dynamic zones which may be updated following mesh motion
    and/or topology change.

See also
    Foam::zoneGeneratorList
    Foam::zoneGenerator
    Foam::zoneSet

SourceFiles
    zonesGenerator.C

\*---------------------------------------------------------------------------*/
#ifndef zonesGenerator_H
#define zonesGenerator_H

#include "zoneGeneratorList.H"
#include "DemandDrivenMeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class zonesGenerator Declaration
\*---------------------------------------------------------------------------*/

class zonesGenerator
:
    public DemandDrivenMeshObject
    <
        polyMesh,
        TopoChangeableMeshObject,
        zonesGenerator
    >,
    public dictionary,
    public zoneGeneratorList
{
    // Private Member Functions

        //- Create IO object for an optional zonesGenerator
        IOobject io(const polyMesh& mesh) const;

        void generate();


protected:

    // Protected Constructors

        friend class DemandDrivenMeshObject;

        //- Read construct from polyMesh
        explicit zonesGenerator(const polyMesh& mesh);


public:

    //- Runtime type information
    TypeName("zonesGenerator");


    // Constructors

        //- Disallow default bitwise copy construction
        zonesGenerator(const zonesGenerator&) = delete;

        using DemandDrivenMeshObject
        <
            polyMesh,
            TopoChangeableMeshObject,
            zonesGenerator
        >::New;


    //- Destructor
    virtual ~zonesGenerator()
    {}


    // Member Functions

        //- Declare zonesGenerator to be a global dictionary
        virtual bool global() const
        {
            return true;
        }

        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        bool writeData(Ostream&) const;

        //- Update for mesh motion
        virtual bool movePoints();

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap& map);

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap& map);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const zonesGenerator&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
