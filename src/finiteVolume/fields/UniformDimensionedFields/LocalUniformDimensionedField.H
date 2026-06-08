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
    Foam::LocalUniformDimensionedField

Description
    Local variant of UniformDimensionedField. Identical to
    UniformDimensionedField, apart from the fact that IO is processor local.

SourceFiles
    LocalUniformDimensionedField.C

\*---------------------------------------------------------------------------*/

#ifndef LocalUniformDimensionedField_H
#define LocalUniformDimensionedField_H

#include "UniformDimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class UniformDimensionedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LocalUniformDimensionedField
:
    public UniformDimensionedField<Type>
{
public:

    //- Runtime type information
    TypeName("LocalUniformDimensionedField");


    // Constructors

        //- Construct from IOobject
        LocalUniformDimensionedField(const IOobject&, const bool read = true);

        //- Construct from IOobject. Either reads or uses supplied value.
        LocalUniformDimensionedField
        (
            const IOobject&,
            const dimensioned<Type>&,
            const bool read = true
        );

        //- Inherit other constructors
        using UniformDimensionedField<Type>::UniformDimensionedField;


    //- Destructor
    virtual ~LocalUniformDimensionedField();


    // Member Functions

        //- Is object global
        virtual bool global() const
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LocalUniformDimensionedField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
