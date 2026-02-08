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
    Foam::uHomogeneousbHomogeneous

Description
    Homogeneous unburnt gas to homogeneous burnt gas composition mapping

SourceFiles
    uHomogeneousbHomogeneous.C

\*---------------------------------------------------------------------------*/

#ifndef uHomogeneousbHomogeneous_H
#define uHomogeneousbHomogeneous_H

#include "ubMixtureMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ubMixtureMaps
{

/*---------------------------------------------------------------------------*\
                   Class uHomogeneousbHomogeneous Declaration
\*---------------------------------------------------------------------------*/

class uHomogeneousbHomogeneous
:
    public ubMixtureMap
{

public:

    //- Runtime type information
    TypeName("uHomogeneousbHomogeneous");


    // Constructors

        //- Construct from components
        uHomogeneousbHomogeneous
        (
            const uRhoMulticomponentThermo& uThermo,
            const bRhoMulticomponentThermo& bThermo
        );

        //- Disallow default bitwise copy construction
        uHomogeneousbHomogeneous(const uHomogeneousbHomogeneous&) = delete;


    //- Destructor
    virtual ~uHomogeneousbHomogeneous();


    // Member Functions

        virtual PtrList<volScalarField::Internal> prompt
        (
            const PtrList<volScalarField>& Yu
        ) const
        {
            NotImplemented;
            return PtrList<volScalarField::Internal>();
        }

        virtual void reset
        (
            const volScalarField& b,
            UPtrList<volScalarField>& Yu,
            const volScalarField& c,
            const UPtrList<const volScalarField>& Yb
        ) const
        {
            FatalErrorInFunction
                << "EGR not supported by " << typeName
                << exit(FatalError);
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const uHomogeneousbHomogeneous&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ubMixtureMaps
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
