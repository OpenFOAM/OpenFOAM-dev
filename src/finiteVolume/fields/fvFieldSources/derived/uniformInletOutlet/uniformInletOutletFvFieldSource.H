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
    Foam::uniformInletOutletFvFieldSource

Description
    This source condition provides a uniform fixed value when the source is
    positive, and the internal value when it is negative (i.e., a sink)

Usage
    \table
        Property          | Description         | Required | Default value
        uniformInletValue | uniform inlet value | yes      |
    \endtable

    Example of the source condition specification with a constant value:
    \verbatim
    <patchName>
    {
        type            uniformInletOutlet;
        uniformInletValue 0.1;
    }
    \endverbatim

    Example of the source condition specification with a time-varying value:
    \verbatim
    <patchName>
    {
        type            uniformInletOutlet;
        uniformInletValue
        {
            type            table;
            values
            (
                (0 0)
                (1 0.1)
                (9 0.1)
                (10 0)
            );
        }
    }
    \endverbatim

SourceFiles
    uniformInletOutletFvFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef uniformInletOutletFvFieldSource_H
#define uniformInletOutletFvFieldSource_H

#include "fvFieldSource.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class uniformInletOutletFvFieldSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class uniformInletOutletFvFieldSource
:
    public fvFieldSource<Type>
{
private:

    // Private Data

        //- Uniform inlet value
        autoPtr<Function1<Type>> uniformInletValue_;


public:

    //- Runtime type information
    TypeName("uniformInletOutlet");


    // Constructors

        //- Construct from internal field and dictionary
        uniformInletOutletFvFieldSource
        (
            const DimensionedField<Type, volMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        uniformInletOutletFvFieldSource
        (
            const uniformInletOutletFvFieldSource<Type>&,
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
                new uniformInletOutletFvFieldSource<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~uniformInletOutletFvFieldSource();


    // Member Functions

        //- Return the source value
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
    #include "uniformInletOutletFvFieldSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
