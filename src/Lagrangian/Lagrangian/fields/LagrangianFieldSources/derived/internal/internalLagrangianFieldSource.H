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
    Foam::internalLagrangianFieldSource

Description
    This source condition retains the internal value.

Usage
    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type            internal;
    }
    \endverbatim

SourceFiles
    internalLagrangianFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef internalLagrangianFieldSource_H
#define internalLagrangianFieldSource_H

#include "LagrangianFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class internalLagrangianFieldSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class internalLagrangianFieldSource
:
    public LagrangianFieldSource<Type>
{
public:

    //- Runtime type information
    TypeName("internal");


    // Constructors

        //- Inherit constructors
        using LagrangianFieldSource<Type>::LagrangianFieldSource;

        //- Construct and return a clone setting the internal field reference
        virtual autoPtr<LagrangianFieldSource<Type>> clone
        (
            const regIOobject& iIo
        ) const
        {
            return autoPtr<LagrangianFieldSource<Type>>
            (
                new internalLagrangianFieldSource<Type>(*this, iIo)
            );
        }


    //- Destructor
    virtual ~internalLagrangianFieldSource();


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "internalLagrangianFieldSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
