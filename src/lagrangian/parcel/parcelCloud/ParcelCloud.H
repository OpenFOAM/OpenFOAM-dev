/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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
    Foam::ParcelCloud

Description
    Outermost template for parcel clouds. Adds the parcelCloud virtualisation
    layer and forwards the methods required by that layer.

\*---------------------------------------------------------------------------*/

#ifndef ParcelCloud_H
#define ParcelCloud_H

#include "parcelCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
class ParcelCloud
:
    public CloudType,
    virtual public parcelCloud
{
public:

    // Constructors

        //- Inherit base constructors
        using CloudType::CloudType;


    // Member Functions

        // Mesh changes

            //- Call this before a topology change. Stores the particles global
            //  positions in the database for use during mapping.
            virtual void storeGlobalPositions()
            {
                CloudType::storeGlobalPositions();
            }

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap& map)
            {
                CloudType::topoChange(map);
            }

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap& map)
            {
                CloudType::mapMesh(map);
            }

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap& map)
            {
                CloudType::distribute(map);
            }


        // Evolution

            //- Evolve the cloud
            virtual void evolve()
            {
                CloudType::evolve();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
