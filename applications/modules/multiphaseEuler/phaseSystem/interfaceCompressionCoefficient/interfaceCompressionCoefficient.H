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
    Foam::interfaceCompressionCoefficient

Description
    ...

SourceFiles
    interfacialCoefficient.C
    interfacialCoefficientNew.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceCompressionCoefficient_H
#define interfaceCompressionCoefficient_H

#include "phaseInterface.H"
#include "BlendedInterfacialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class interfaceCompressionCoefficient Declaration
\*---------------------------------------------------------------------------*/

class interfaceCompressionCoefficient
:
    public dimensionedScalar
{
public:

    //- Runtime type information
    TypeName("interfaceCompressionCoefficient");


    // Constructors

        //- Inherit constructors
        using dimensionedScalar::dimensionedScalar;

        //- Disallow default bitwise copy construction
        interfaceCompressionCoefficient
        (
            const interfaceCompressionCoefficient&
        ) = delete;


    //- Destructor
    virtual ~interfaceCompressionCoefficient();


    // Selectors

        static autoPtr<interfaceCompressionCoefficient> New
        (
            const UPtrList<const entry>& entries,
            const phaseInterface& interface
        );


    // Member Functions

        //- Return the interface compression coefficient
        const dimensionedScalar& cAlpha() const;
};


/*---------------------------------------------------------------------------*\
            Class blendedInterfaceCompressionCoefficient Declaration
\*---------------------------------------------------------------------------*/

class blendedInterfaceCompressionCoefficient
:
    public BlendedInterfacialModel<interfaceCompressionCoefficient>
{
public:

    // Constructors

        //- Inherit base class constructors
        using
            BlendedInterfacialModel<interfaceCompressionCoefficient>::
            BlendedInterfacialModel;


    // Selectors

        inline static autoPtr<blendedInterfaceCompressionCoefficient> New
        (
            const UPtrList<const entry>& entries,
            const phaseInterface& interface,
            const entry& blendingEntry
        )
        {
            return autoPtr<blendedInterfaceCompressionCoefficient>
            (
                new blendedInterfaceCompressionCoefficient
                (
                    entries,
                    interface,
                    blendingEntry
                )
            );
        }


    // Member Functions

        //- Return the interface compression coefficient
        template<class GeoMesh>
        tmp<GeometricField<scalar, GeoMesh>> cAlpha() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
