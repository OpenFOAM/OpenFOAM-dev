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
    Foam::carrierLagrangianFieldSource

Description
    This source condition creates a Lagrangian value interpolated from a given
    carrier field. By default the field chosen is that with the same name as
    the Lagrangian field, but this can be overridden with the \c fieldc
    control.

Usage
    \table
        Property  | Description            | Required? | Default
        fieldc    | The carrier field name | no | The Lagrangian field name
    \endtable

    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type            carrier;
        fieldc          U.air;
    }
    \endverbatim

SourceFiles
    carrierLagrangianFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef carrierLagrangianFieldSource_H
#define carrierLagrangianFieldSource_H

#include "LagrangianFieldSource.H"
#include "CloudLagrangianFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class carrierLagrangianFieldSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class carrierLagrangianFieldSource
:
    public LagrangianFieldSource<Type>,
    private CloudLagrangianFieldSource<Type>
{
private:

    // Private Data

        //- The carrier field name
        const word fieldcName_;


public:

    //- Runtime type information
    TypeName("carrier");


    // Constructors

        //- Construct from internal field and dictionary
        carrierLagrangianFieldSource
        (
            const regIOobject&,
            const dictionary& dict
        );

        //- Copy constructor setting the internal field reference
        carrierLagrangianFieldSource
        (
            const carrierLagrangianFieldSource&,
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
                new carrierLagrangianFieldSource(*this, iIo)
            );
        }


    //- Destructor
    virtual ~carrierLagrangianFieldSource();


    // Member Functions

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
    #include "carrierLagrangianFieldSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
