/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::parcelClouds

Description
    List of parcel clouds, with the same interface as an individual parcel
    cloud. Is a mesh object, so mesh change hooks are provided and will be
    applied to the contained cloud. This is the object that should be
    constructed by a solver in order to support the coupled simulation of
    multiple clouds. An fvModel should *not* construct this object, as that
    would nest two mesh objects. An fvModel should construct the base
    parcelCloudList instead.

SourceFiles
    parcelClouds.C

\*---------------------------------------------------------------------------*/

#ifndef parcelClouds_H
#define parcelClouds_H

#include "DemandDrivenMeshObject.H"
#include "parcelCloudList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                          Class parcelClouds Declaration
\*---------------------------------------------------------------------------*/

class parcelClouds
:
    public DemandDrivenMeshObject<fvMesh, UpdateableMeshObject, parcelClouds>,
    public parcelCloudList
{
protected:

    friend class DemandDrivenMeshObject
    <
        fvMesh,
        UpdateableMeshObject,
        parcelClouds
    >;

    // Protected Constructors

        //- Construct with given mesh and carrier fields
        parcelClouds
        (
            const fvMesh& mesh,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g
        );

        //- Construct with given mesh and carrier fields and thermo
        parcelClouds
        (
            const fvMesh& mesh,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const fluidThermo& carrierThermo
        );


public:

    // Constructors

        //- Disallow default bitwise copy construction
        parcelClouds(const parcelClouds&) = delete;

        //- Inherit the base New method
        using DemandDrivenMeshObject
        <
            fvMesh,
            UpdateableMeshObject,
            parcelClouds
        >::New;


    //- Destructor
    virtual ~parcelClouds();


    // Member Functions

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
