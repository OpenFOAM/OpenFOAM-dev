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
    Foam::uniformFixedValueLagrangianFieldSource

Description
    This source condition provides a uniform fixed value.

Usage
    \table
        Property          | Description         | Required? | Default
        uniformValue      | The uniform value   | yes       |
    \endtable

    Example specification for a constant value:
    \verbatim
    <LagrangianModelName>
    {
        type            uniformFixedValue;
        uniformValue    0.1;
    }
    \endverbatim

    Example specification for a time-varying value:
    \verbatim
    <LagrangianModelName>
    {
        type            uniformFixedValue;
        uniformValue
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
    uniformFixedValueLagrangianFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef uniformFixedValueLagrangianFieldSource_H
#define uniformFixedValueLagrangianFieldSource_H

#include "LagrangianFieldSource.H"
#include "Function1LagrangianFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class uniformFixedValueLagrangianFieldSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class uniformFixedValueLagrangianFieldSource
:
    public LagrangianFieldSource<Type>,
    private Function1LagrangianFieldSource<Type>
{
private:

    // Private Data

        //- Uniform value
        autoPtr<Function1<Type>> uniformValue_;


    // Private Member Functions

        //- Use the base class' value method
        using Function1LagrangianFieldSource<Type>::value;


public:

    //- Runtime type information
    TypeName("uniformFixedValue");


    // Constructors

        //- Construct from internal field and dictionary
        uniformFixedValueLagrangianFieldSource
        (
            const regIOobject&,
            const dictionary& dict
        );

        //- Copy constructor setting the internal field reference
        uniformFixedValueLagrangianFieldSource
        (
            const uniformFixedValueLagrangianFieldSource<Type>&,
            const regIOobject&
        );

        //- Construct and return a clone setting the internal field reference
        virtual autoPtr<LagrangianFieldSource<Type>> clone
        (
            const regIOobject& iIo
        ) const
        {
            return autoPtr<LagrangianFieldSource<Type>>
            (
                new uniformFixedValueLagrangianFieldSource<Type>(*this, iIo)
            );
        }


    //- Destructor
    virtual ~uniformFixedValueLagrangianFieldSource();


    // Member Functions

        //- Return the source value
        virtual tmp<LagrangianSubField<Type>> sourceValue
        (
            const LagrangianSource&,
            const LagrangianSubMesh&
        ) const;

        //- Return the internal coefficient
        virtual tmp<LagrangianSubScalarField> internalCoeff
        (
            const LagrangianSource&,
            const LagrangianSubMesh&
        ) const;

        //- Return the value for an instantaneous injection
        virtual tmp<LagrangianSubField<Type>> value
        (
            const LagrangianInjection&,
            const LagrangianSubMesh&
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "uniformFixedValueLagrangianFieldSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
