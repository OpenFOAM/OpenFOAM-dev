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
    Foam::genericLagrangianFieldSource

Description
    This provides a generic source condition, useful as a fallback for handling
    unknown types when post-processing or running mesh manipulation utilities.
    Not generally applicable as a user-specified condition.

SourceFiles
    genericLagrangianFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef genericLagrangianFieldSource_H
#define genericLagrangianFieldSource_H

#include "genericFieldBase.H"
#include "LagrangianFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class genericLagrangianFieldSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class genericLagrangianFieldSource
:
    public genericFieldBase,
    public LagrangianFieldSource<Type>
{
    // Private Data

        //- Field source dictionary
        dictionary dict_;


public:

    //- Runtime type information
    TypeName("generic");


    // Constructors

        //- Construct from internal field and dictionary
        genericLagrangianFieldSource
        (
            const regIOobject&,
            const dictionary&
        );

        //- Construct as copy setting internal field reference
        genericLagrangianFieldSource
        (
            const genericLagrangianFieldSource<Type>&,
            const regIOobject&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<LagrangianFieldSource<Type>> clone
        (
            const regIOobject& iIo
        ) const
        {
            return autoPtr<LagrangianFieldSource<Type>>
            (
                new genericLagrangianFieldSource<Type>(*this, iIo)
            );
        }


    //- Destructor
    virtual ~genericLagrangianFieldSource();


    // Member Functions

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "genericLagrangianFieldSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
