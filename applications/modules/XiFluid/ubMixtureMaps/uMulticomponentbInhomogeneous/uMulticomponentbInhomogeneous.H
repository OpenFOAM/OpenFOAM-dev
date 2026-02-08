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
    Foam::uMulticomponentbInhomogeneous

Description
    Multicomponent unburnt gas to inhomogeneous burnt gas composition mapping

SourceFiles
    uMulticomponentbInhomogeneous.C

\*---------------------------------------------------------------------------*/

#ifndef uMulticomponentbInhomogeneous_H
#define uMulticomponentbInhomogeneous_H

#include "ubMixtureMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ubMixtureMaps
{

/*---------------------------------------------------------------------------*\
                   Class uMulticomponentbInhomogeneous Declaration
\*---------------------------------------------------------------------------*/

class uMulticomponentbInhomogeneous
:
    public ubMixtureMap
{
    // Private Data

        //- Fuel specie index
        label fu_;

        //- Oxidant species mass fractions
        List<Tuple2<label, scalar>> ox_;

        //- Product species mass fractions
        List<Tuple2<label, scalar>> pr_;


public:

    //- Runtime type information
    TypeName("uMulticomponentbInhomogeneous");


    // Constructors

        //- Construct from components
        uMulticomponentbInhomogeneous
        (
            const uRhoMulticomponentThermo& uThermo,
            const bRhoMulticomponentThermo& bThermo
        );

        //- Disallow default bitwise copy construction
        uMulticomponentbInhomogeneous
        (
            const uMulticomponentbInhomogeneous&
        ) = delete;


    //- Destructor
    virtual ~uMulticomponentbInhomogeneous();


    // Member Functions

        virtual PtrList<volScalarField::Internal> prompt
        (
            const PtrList<volScalarField>& Yu
        ) const;

        virtual void reset
        (
            const volScalarField& b,
            UPtrList<volScalarField>& Yu,
            const volScalarField& c,
            const UPtrList<const volScalarField>& Yb
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const uMulticomponentbInhomogeneous&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ubMixtureMaps
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
