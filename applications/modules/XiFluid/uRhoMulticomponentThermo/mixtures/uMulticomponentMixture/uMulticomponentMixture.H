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
    Foam::uMulticomponentMixture

Description
    Unburnt gas combustion mixture for inhomogeneous lean, stoichiometric
    or rich conditions

    in which the fuel mass-fraction \c fu defines the composition.

SourceFiles
    uMulticomponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef uMulticomponentMixture_H
#define uMulticomponentMixture_H

#include "coefficientMulticomponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class uMulticomponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class uMulticomponentMixture
:
    public coefficientMulticomponentMixture<ThermoType>
{

public:

    // Constructors

        //- Construct from a dictionary
        uMulticomponentMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        uMulticomponentMixture
        (
            const uMulticomponentMixture<ThermoType>&
        ) = delete;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "uMulticomponentMixture<" + ThermoType::typeName() + '>';
        }

        //- Return the fuel-oxidant equivalence ratio
        scalar Phi(const scalarFieldListSlice& Y) const;

        //- Return the prompt burnt gas composition
        PtrList<volScalarField::Internal> prompt
        (
            const PtrList<volScalarField>& Yu
        ) const;

        void reset
        (
            const volScalarField& b,
            PtrList<volScalarField>& Yu,
            const volScalarField& c,
            const PtrList<volScalarField>& Yb
        ) const
        {
            NotImplemented;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const uMulticomponentMixture<ThermoType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "uMulticomponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
