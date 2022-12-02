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

        //- The vol field to interpolate
        const VolField<Type>& psi_;

        //- Reference to the mesh
        const polyMesh& mesh_;


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
                const VolField<Type>& psi
            ),
            (psi)
        );


    // Selectors

        //- Return a reference to the specified interpolation scheme
        static autoPtr<interpolation<Type>> New
        (
            const word& interpolationType,
            const VolField<Type>& psi
        );

        //- Return a reference to the selected interpolation scheme
        static autoPtr<interpolation<Type>> New
        (
            const dictionary& interpolationSchemes,
            const VolField<Type>& psi
        );


    // Constructors

        //- Construct from components
        interpolation
        (
            const VolField<Type>& psi
        );


    //- Destructor
    virtual ~interpolation()
    {}


    // Member Functions

        //- Return the field to be interpolated
        const VolField<Type>& psi() const
        {
            return psi_;
        }

        //- Interpolate field to the given point in the given cell
        virtual Type interpolate
        (
            const vector& position,
            const label celli,
            const label facei = -1
        ) const = 0;

        //- As above, but for a field
        virtual tmp<Field<Type>> interpolate
        (
            const vectorField& position,
            const labelList& celli,
            const labelList& facei = NullObjectRef<labelList>()
        ) const = 0;

        //- Interpolate field to the given coordinates in the tetrahedron
        //  defined by the given indices. Calls interpolate function
        //  above here except where overridden by derived
        //  interpolation types.
        virtual Type interpolate
        (
            const barycentric& coordinates,
            const tetIndices& tetIs,
            const label facei = -1
        ) const;

        //- As above, but for a field
        virtual tmp<Field<Type>> interpolate
        (
            const Field<barycentric>& coordinates,
            const labelList& celli,
            const labelList& tetFacei,
            const labelList& tetPti,
            const labelList& facei = NullObjectRef<labelList>()
        ) const = 0;
};


/*---------------------------------------------------------------------------*\
                     Class fieldInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class InterpolationType>
class fieldInterpolation
:
    public interpolation<Type>
{
public:

    // Constructors

        //- Inherit base class constructors
        using interpolation<Type>::interpolation;


    // Member Functions

        //- Interpolate field to the given points in the given cells
        virtual tmp<Field<Type>> interpolate
        (
            const vectorField& position,
            const labelList& celli,
            const labelList& facei = NullObjectRef<labelList>()
        ) const;

        //- Interpolate field to the given coordinates in the given tetrahedra
        virtual tmp<Field<Type>> interpolate
        (
            const Field<barycentric>& coordinates,
            const labelList& celli,
            const labelList& tetFacei,
            const labelList& tetPti,
            const labelList& facei = NullObjectRef<labelList>()
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineInterpolation(Type, InterpolationType)                           \
    defineNamedTemplateTypeNameAndDebug(InterpolationType<Type>, 0);

#define makeInterpolation(Type, InterpolationType)                             \
    interpolation<Type>::                                                      \
        adddictionaryConstructorToTable<InterpolationType<Type>>               \
        add##InterpolationType##Type##ConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
