/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::volPointInterpolation

Description
    Interpolate from cell centres to points (vertices) using inverse distance
    weighting

SourceFiles
    volPointInterpolation.C
    volPointInterpolationTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef volPointInterpolation_H
#define volPointInterpolation_H

#include "DemandDrivenMeshObject.H"
#include "scalarList.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class volPointInterpolation Declaration
\*---------------------------------------------------------------------------*/

class volPointInterpolation
:
    public DemandDrivenMeshObject
    <
        fvMesh,
        UpdateableMeshObject,
        volPointInterpolation
    >
{
    // Private Data

        //- Interpolation scheme weighting factor array.
        scalarListList pointWeights_;

        //- Boundary addressing
        autoPtr<primitivePatch> boundaryPtr_;

        //- Per boundary point the weights per pointFaces.
        List<scalarListList> boundaryPointWeights_;

        //- Per boundary point the weights per pointFaces.
        List<scalarListList> boundaryPointNbrWeights_;


    // Private Member Functions

        //- Construct all point weighting factors
        void makeWeights();


protected:

    friend class DemandDrivenMeshObject
    <
        fvMesh,
        UpdateableMeshObject,
        volPointInterpolation
    >;

    // Protected Constructors

        //- Constructor given fvMesh
        explicit volPointInterpolation(const fvMesh&);


public:

    // Declare name of the class and its debug switch
    ClassName("volPointInterpolation");


    // Constructors

        //- Disallow default bitwise copy construction
        volPointInterpolation(const volPointInterpolation&) = delete;


    //- Destructor
    ~volPointInterpolation();


    // Member Functions

        //- Correct weighting factors for moving mesh.
        virtual bool movePoints();

        //- Update mesh topology using the morph engine
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Interpolation functions

        //- Interpolate volField using inverse distance weighting
        //  returning pointField
        template<class Type>
        tmp<PointField<Type>> interpolate
        (
            const VolField<Type>&
        ) const;

        //- Interpolate tmp<volField> using inverse distance weighting
        //  returning pointField
        template<class Type>
        tmp<PointField<Type>> interpolate
        (
            const tmp<VolField<Type>>&
        ) const;


        // Low level

            //- Interpolate from volField to pointField
            //  using inverse distance weighting
            template<class Type>
            void interpolateUnconstrained
            (
                const VolField<Type>&,
                PointField<Type>&
            ) const;

            //- Interpolate from volField to pointField
            //  using inverse distance weighting
            template<class Type>
            void interpolate
            (
                const VolField<Type>&,
                PointField<Type>&
            ) const;

            //- Interpolate volField using inverse distance weighting
            //  returning pointField with name. Optionally caches.
            template<class Type>
            tmp<PointField<Type>> interpolate
            (
                const VolField<Type>&,
                const word& name,
                const bool cache
            ) const;


        // Interpolation for displacement (applies 2D correction)

            //- Interpolate from volField to pointField
            //  using inverse distance weighting
            void interpolateDisplacement
            (
                const volVectorField&,
                pointVectorField&
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const volPointInterpolation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "volPointInterpolationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
