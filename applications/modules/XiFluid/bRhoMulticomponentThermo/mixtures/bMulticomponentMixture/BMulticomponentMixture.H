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
    Foam::BMulticomponentMixture

Description
    Templated burnt gas combustion multicomponent mixture

See also
    Foam::bMulticomponentMixture

SourceFiles
    BMulticomponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef BMulticomponentMixture_H
#define BMulticomponentMixture_H

#include "coefficientMulticomponentMixture.H"
#include "bMulticomponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class BMulticomponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class BMulticomponentMixture
:
    public coefficientMulticomponentMixture<ThermoType>,
    public bMulticomponentMixture
{

public:

    // Constructors

        //- Construct from a dictionary
        BMulticomponentMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        BMulticomponentMixture
        (
            const BMulticomponentMixture<ThermoType>&
        ) = delete;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                bMulticomponentMixture::typeName
              + '<' + ThermoType::typeName() + '>';
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const BMulticomponentMixture<ThermoType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BMulticomponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
