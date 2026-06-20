/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2026 OpenFOAM Foundation
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
    Foam::MRFZones

Description
    MRF zones DemandDrivenMeshObject read from MRFProperties dictionary

    MRF zones are specified by a list of dictionary entries, e.g.

    \verbatim
    zone1
    {
        cellZone    rotor1;
        ...
    }

    zone2
    {
        cellZone    rotor2;
        ...
    }
    \endverbatim

SourceFiles
    MRFZones.C

\*---------------------------------------------------------------------------*/

#ifndef MRFZones_H
#define MRFZones_H

#include "MRFZoneList.H"
#include "DemandDrivenMeshObject.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class MRFZones Declaration
\*---------------------------------------------------------------------------*/

class MRFZones
:
    public DemandDrivenMeshObject
    <
        fvMesh,
        TopoChangeableMeshObject,
        MRFZones,
        IOdictionary
    >,
    public MRFZoneList
{
    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;


protected:

    friend class DemandDrivenMeshObject
    <
        fvMesh,
        TopoChangeableMeshObject,
        MRFZones,
        IOdictionary
    >;

    // Protected Constructors

        //- Construct from mesh
        explicit MRFZones(const fvMesh& mesh);


public:

    //- Runtime type information
    TypeName("MRFZones");


    // Constructors

        //- Disallow default bitwise copy construction
        MRFZones(const MRFZones&) = delete;

        //- Inherit the base New method
        using DemandDrivenMeshObject
        <
            fvMesh,
            TopoChangeableMeshObject,
            MRFZones,
            IOdictionary
        >::New;


        //- Destructor
        virtual ~MRFZones()
        {}


    // Member Functions

        //- Use MRFZoneList::size for looping
        using MRFZoneList::size;

        // Mesh changes

            //- Prepare for mesh update
            virtual void preUpdateMesh();

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // IO

            //- Read dictionary
            virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const MRFZones&) = delete;
};


/*---------------------------------------------------------------------------*\
                           Class NullMRF Declaration
\*---------------------------------------------------------------------------*/

class NullMRF
{
public:

    // Constructors

        //- Construct null
        NullMRF()
        {}


    // Member Functions

        //- Return the argument unchanged
        template<class Type>
        inline const Type& relative(const Type& U) const
        {
            return U;
        }

        //- Return the argument unchanged
        template<class Type>
        inline const Type& relative(const Type& U, const label patchi) const
        {
            return U;
        }

        //- Return the argument unchanged
        template<class Type>
        inline const Type& absolute(const Type& U) const
        {
            return U;
        }

        //- Return the argument unchanged
        template<class Type>
        inline const Type& absolute(const Type& U, const label patchi) const
        {
            return U;
        }


        //- Return the argument unchanged
        template<class Type, class RhoType>
        inline const Type& relative(const Type& U, const RhoType& rho) const
        {
            return U;
        }

        //- Return the argument unchanged
        template<class Type, class RhoType>
        inline const Type& relative
        (
            const Type& U,
            const RhoType& rho,
            const label patchi
        ) const
        {
            return U;
        }

        //- Return the argument unchanged
        template<class Type, class RhoType>
        inline const Type& absolute(const Type& U, const RhoType& rho) const
        {
            return U;
        }

        //- Return the argument unchanged
        template<class Type, class RhoType>
        inline const Type& absolute
        (
            const Type& U,
            const RhoType& rho,
            const label patchi
        ) const
        {
            return U;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
