/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::NaNFvFieldSource

Description
    This source condition provides a NaN value.

Usage
    Example specification:
    \verbatim
    <sourceName>
    {
        type            NaN;
    }
    \endverbatim

SourceFiles
    NaNFvFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef NaNFvFieldSource_H
#define NaNFvFieldSource_H

#include "fvFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class NaNFvFieldSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class NaNFvFieldSource
:
    public fvFieldSource<Type>
{
public:

    //- Runtime type information
    TypeName("NaN");


    // Constructors

        //- Construct from internal field and dictionary
        NaNFvFieldSource
        (
            const DimensionedField<Type, volMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        NaNFvFieldSource
        (
            const NaNFvFieldSource<Type>&,
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
                new NaNFvFieldSource<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~NaNFvFieldSource();


    // Member Functions

        //- Return the source value
        virtual tmp<DimensionedField<Type, volMesh>> sourceValue
        (
            const fvSource& model,
            const DimensionedField<scalar, volMesh>& source
        ) const;

        //- Return the source value
        virtual tmp<Field<Type>> sourceValue
        (
            const fvSource& model,
            const scalarField& source,
            const labelUList& cells
        ) const;

        //- Return the source value
        virtual tmp<DimensionedField<scalar, volMesh>> internalCoeff
        (
            const fvSource& model,
            const DimensionedField<scalar, volMesh>& source
        ) const;

        //- Return the internal coefficient
        virtual tmp<scalarField> internalCoeff
        (
            const fvSource& model,
            const scalarField& source,
            const labelUList& cells
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NaNFvFieldSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
