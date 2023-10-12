/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::internalFvFieldSource

Description
    This source condition provides the internal value.

Usage
    Example of the source condition specification:
    \verbatim
    <patchName>
    {
        type            internal;
    }
    \endverbatim

SourceFiles
    internalFvFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef internalFvFieldSource_H
#define internalFvFieldSource_H

#include "fvFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class internalFvFieldSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class internalFvFieldSource
:
    public fvFieldSource<Type>
{
public:

    //- Runtime type information
    TypeName("internal");


    // Constructors

        //- Construct from internal field and dictionary
        internalFvFieldSource
        (
            const DimensionedField<Type, volMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        internalFvFieldSource
        (
            const internalFvFieldSource<Type>&,
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
                new internalFvFieldSource<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~internalFvFieldSource();


    // Member Functions

        //- Return the source value
        virtual tmp<Field<Type>> sourceValue(const fvSource&) const;

        //- Return the internal coefficient
        virtual tmp<scalarField> internalCoeff(const fvSource&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "internalFvFieldSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
