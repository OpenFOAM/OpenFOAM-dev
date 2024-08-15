/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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
    Foam::genericFvFieldSource

Description
    This provides a generic source condition, useful as a fallback for handling
    unknown types when post-processing or running mesh manipulation utilities.
    Not generally applicable as a user-specified condition.

SourceFiles
    genericFvFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef genericFvFieldSource_H
#define genericFvFieldSource_H

#include "genericFieldBase.H"
#include "fvFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class genericFvFieldSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class genericFvFieldSource
:
    public genericFieldBase,
    public fvFieldSource<Type>
{
    // Private Data

        //- Field source dictionary
        dictionary dict_;


public:

    //- Runtime type information
    TypeName("generic");


    // Constructors

        //- Construct internal field and dictionary
        genericFvFieldSource
        (
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Copy constructor setting internal field reference
        genericFvFieldSource
        (
            const genericFvFieldSource<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<fvFieldSource<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return autoPtr<fvFieldSource<Type>>
            (
                new genericFvFieldSource<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~genericFvFieldSource();


    // Member Functions

        //- Return the source energy value
        virtual tmp<Field<Type>> sourceValue(const fvSource&) const;

        //- Return the internal coefficient
        virtual tmp<scalarField> internalCoeff(const fvSource&) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "genericFvFieldSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
