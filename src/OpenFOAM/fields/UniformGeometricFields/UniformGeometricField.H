/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::UniformGeometricField

Description
    Dimensioned<Type> registered with the database as a registered IOobject
    which has the functionality of a uniform field and allows values from the
    top-level code to be passed to boundary conditions etc.

SourceFiles
    UniformGeometricField.C

\*---------------------------------------------------------------------------*/

#ifndef UniformGeometricField_H
#define UniformGeometricField_H

#include "UniformDimensionedField.H"
#include "UniformFieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class UniformGeometricField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class UniformGeometricField
:
    public UniformDimensionedField<Type>
{
public:

    // Public Typedefs

        typedef UniformDimensionedField<Type> Internal;
        typedef UniformField<Type> Patch;
        typedef UniformFieldField<Type> Boundary;
        typedef Type cmptType;


    //- Runtime type information
    TypeName("UniformGeometricField");


    // Constructors

        //- Construct from components. Either reads or uses supplied value.
        UniformGeometricField(const IOobject&, const dimensioned<Type>&);

        //- Copy constructor
        UniformGeometricField(const UniformGeometricField<Type>&);

        //- Construct from Istream
        UniformGeometricField(const IOobject&);


    //- Destructor
    virtual ~UniformGeometricField();


    // Member Operators

        inline const UniformGeometricField<Type>& oldTime() const;

        inline const Internal& operator()() const;

        inline const Internal& v() const;

        inline typename Internal::FieldType primitiveField() const;

        inline Boundary boundaryField() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UniformGeometricField.C"
    #include "UniformGeometricFieldI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
