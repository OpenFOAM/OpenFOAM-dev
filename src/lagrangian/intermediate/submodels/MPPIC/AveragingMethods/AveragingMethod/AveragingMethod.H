/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::AveragingMethod

Description
    Base class for lagrangian averaging methods.

SourceFiles
    AveragingMethod.C
    AveragingMethodI.H

\*---------------------------------------------------------------------------*/

#ifndef AveragingMethod_H
#define AveragingMethod_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "barycentric.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class AveragingMethod Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class AveragingMethod
:
    public regIOobject,
    public FieldField<Field, Type>
{
protected:

    //- Protected typedefs

        //- Gradient type
        typedef typename outerProduct<vector, Type>::type TypeGrad;


    //- Protected data

        //- Dictionary
        const dictionary& dict_;

        //- The mesh on which the averaging is to be done
        const fvMesh& mesh_;


    //- Protected member functions

        //- Update the gradient calculation
        virtual void updateGrad();


public:

    //- Runtime type information
    TypeName("averagingMethod");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        AveragingMethod,
        dictionary,
        (
            const IOobject& io,
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (io, dict, mesh)
    );


    //- Constructors

        //- Construct from components
        AveragingMethod
        (
            const IOobject& io,
            const dictionary& dict,
            const fvMesh& mesh,
            const labelList& size
        );

        //- Construct a copy
        AveragingMethod(const AveragingMethod<Type>& am);

        //- Construct and return a clone
        virtual autoPtr<AveragingMethod<Type>> clone() const = 0;


    //- Selector
    static autoPtr<AveragingMethod<Type>> New
    (
        const IOobject& io,
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~AveragingMethod();


    //- Member Functions

        //- Add point value to interpolation
        virtual void add
        (
            const barycentric& coordinates,
            const tetIndices& tetIs,
            const Type& value
        ) = 0;

        //- Interpolate
        virtual Type interpolate
        (
            const barycentric& coordinates,
            const tetIndices& tetIs
        ) const = 0;

        //- Interpolate gradient
        virtual TypeGrad interpolateGrad
        (
            const barycentric& coordinates,
            const tetIndices& tetIs
        ) const = 0;

        //- Calculate the average
        virtual void average();
        virtual void average(const AveragingMethod<scalar>& weight);

        //- Dummy write
        virtual bool writeData(Ostream&) const;

        //- Write using setting from DB
        virtual bool write(const bool valid = true) const;

        //- Return an internal field of the average
        virtual tmp<Field<Type>> primitiveField() const = 0;

        //- Assign to another average
        inline void operator=(const AveragingMethod<Type>& x);

        //- Assign to value
        inline void operator=(const Type& x);

        //- Assign to tmp
        inline void operator=(tmp<FieldField<Field, Type>> x);

        //- Add-equal tmp
        inline void operator+=(tmp<FieldField<Field, Type>> x);

        //- Multiply-equal tmp
        inline void operator*=(tmp<FieldField<Field, Type>> x);

        //- Divide-equal tmp
        inline void operator/=(tmp<FieldField<Field, scalar>> x);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AveragingMethodI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AveragingMethod.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
