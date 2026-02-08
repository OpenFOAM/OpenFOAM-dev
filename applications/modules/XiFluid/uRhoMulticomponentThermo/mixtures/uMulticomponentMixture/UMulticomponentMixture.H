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
    Foam::UMulticomponentMixture

Description
    Templated unburnt gas multicomponent mixture

See also
    Foam::uMulticomponentMixture

SourceFiles
    UMulticomponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef UMulticomponentMixture_H
#define UMulticomponentMixture_H

#include "coefficientMulticomponentMixture.H"
#include "uMulticomponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class UMulticomponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class UMulticomponentMixture
:
    public coefficientMulticomponentMixture<ThermoType>,
    public uMulticomponentMixture
{

public:

    // Constructors

        //- Construct from a dictionary
        UMulticomponentMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        UMulticomponentMixture
        (
            const UMulticomponentMixture<ThermoType>&
        ) = delete;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                uMulticomponentMixture::typeName
              + '<' + ThermoType::typeName() + '>';
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const UMulticomponentMixture<ThermoType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UMulticomponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
