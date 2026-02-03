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
    Foam::bMulticomponentMixture

Description
    Burnt gas combustion multicomponent mixture

SourceFiles
    bMulticomponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef bMulticomponentMixture_H
#define bMulticomponentMixture_H

#include "coefficientMulticomponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class bMulticomponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class bMulticomponentMixture
:
    public coefficientMulticomponentMixture<ThermoType>
{

public:

    // Constructors

        //- Construct from a dictionary
        bMulticomponentMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        bMulticomponentMixture
        (
            const bMulticomponentMixture<ThermoType>&
        ) = delete;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "bMulticomponentMixture<" + ThermoType::typeName() + '>';
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const bMulticomponentMixture<ThermoType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "bMulticomponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
