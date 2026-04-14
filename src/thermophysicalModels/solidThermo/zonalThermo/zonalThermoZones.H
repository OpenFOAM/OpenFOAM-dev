/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::zonalThermoZones

Description
    Mesh object to store cell-zone correspondence for zonal thermo.

SourceFiles
    zonalThermoZones.C

\*---------------------------------------------------------------------------*/

#ifndef zonalThermoZones_H
#define zonalThermoZones_H

#include "DemandDrivenMeshObject.H"
#include "fvMesh.H"
#include "hashedWordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class zonalThermoZones Declaration
\*---------------------------------------------------------------------------*/

class zonalThermoZones
:
    public DemandDrivenMeshObject
    <
        fvMesh,
        TopoChangeableMeshObject,
        zonalThermoZones
    >
{
    // Private Data

        //- The names of the zones
        const hashedWordList& zones_;

        //- For each cell, the (local) index of the zone
        labelList cellZones_;


    // Private Member Functions

        //- Update the cell-zones map
        void update();


protected:

    friend class DemandDrivenMeshObject
    <
        fvMesh,
        TopoChangeableMeshObject,
        zonalThermoZones
    >;


    // Protected Constructors

        //- Construct from a name, a mesh and a list of zones
        explicit zonalThermoZones
        (
            const word& name,
            const fvMesh& mesh,
            const hashedWordList& zones
        );


public:

    //- Runtime type information
    TypeName("zonalThermoZones");


    // Constructors

        //- Disallow default bitwise copy construction
        zonalThermoZones(const zonalThermoZones&) = delete;

        //- Inherit the base New method
        using DemandDrivenMeshObject
        <
            fvMesh,
            TopoChangeableMeshObject,
            zonalThermoZones
        >::New;


    //- Destructor
    virtual ~zonalThermoZones();


    // Member Functions

        //- Access the cell zone indices
        inline const labelList& cellZones() const
        {
            return cellZones_;
        }

        //- Update for mesh motion
        virtual bool movePoints();

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
