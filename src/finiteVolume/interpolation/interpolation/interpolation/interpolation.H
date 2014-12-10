/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::interpolation

Description
    Abstract base class for interpolation

\*---------------------------------------------------------------------------*/

#ifndef interpolation_H
#define interpolation_H

#include "faceList.H"
#include "volFieldsFwd.H"
#include "pointFields.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "tetIndices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class interpolation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolation
{

protected:

    // Protected data

        const GeometricField<Type, fvPatchField, volMesh>& psi_;

        const polyMesh& pMesh_;
        const vectorField& pMeshPoints_;
        const faceList& pMeshFaces_;
        const vectorField& pMeshFaceCentres_;
        const vectorField& pMeshFaceAreas_;


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            interpolation,
            dictionary,
            (
                const GeometricField<Type, fvPatchField, volMesh>& psi
            ),
            (psi)
        );


    // Selectors

        //- Return a reference to the specified interpolation scheme
        static autoPtr<interpolation<Type> > New
        (
            const word& interpolationType,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        //- Return a reference to the selected interpolation scheme
        static autoPtr<interpolation<Type> > New
        (
            const dictionary& interpolationSchemes,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );


    // Constructors

        //- Construct from components
        interpolation
        (
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );


    //- Destructor
    virtual ~interpolation()
    {}


    // Member Functions

        //- Return the field to be interpolated
        const GeometricField<Type, fvPatchField, volMesh>& psi() const
        {
            return psi_;
        }

        //- Interpolate field to the given point in the given cell
        virtual Type interpolate
        (
            const vector& position,
            const label cellI,
            const label faceI = -1
        ) const = 0;

        //- Interpolate field to the given point in the tetrahedron
        //  defined by the given indices.  Calls interpolate function
        //  above here execpt where overridden by derived
        //  interpolation types.
        virtual Type interpolate
        (
            const vector& position,
            const tetIndices& tetIs,
            const label faceI = -1
        ) const
        {
            return interpolate(position, tetIs.cell(), faceI);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeInterpolationType(SS, Type)                                        \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
interpolation<Type>::adddictionaryConstructorToTable<SS<Type> >                \
    add##SS##Type##ConstructorToTable_;


#define makeInterpolation(SS)                                                  \
                                                                               \
makeInterpolationType(SS, scalar)                                              \
makeInterpolationType(SS, vector)                                              \
makeInterpolationType(SS, sphericalTensor)                                     \
makeInterpolationType(SS, symmTensor)                                          \
makeInterpolationType(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "interpolation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
