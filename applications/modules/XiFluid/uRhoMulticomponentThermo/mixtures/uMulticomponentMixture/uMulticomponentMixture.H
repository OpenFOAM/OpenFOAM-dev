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
    Unburnt gas combustion multicomponent mixture

SourceFiles
    uMulticomponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef uMulticomponentMixture_H
#define uMulticomponentMixture_H

#include "uMixture.H"
#include "speciesTable.H"
#include "FieldListSlice.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class uMulticomponentMixture Declaration
\*---------------------------------------------------------------------------*/

class uMulticomponentMixture
:
    public uMixture
{
    // Private Data

        //- Fuel specie index
        label fu_;

        //- Stoichiometric air/fuel ratio
        scalar stoicRatio_;


public:

    //- Runtime type information
    ClassName("uMulticomponent");


    // Constructors

        //- Construct from the speciesTable and dictionary
        uMulticomponentMixture
        (
            const speciesTable& species,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        uMulticomponentMixture(const uMulticomponentMixture&) = delete;


    // Member Functions

        //- Return the mixture type name
        virtual const word& mixtureType() const
        {
            return typeName;
        }

        //- Return the fuel specie index
        label fu() const
        {
            return fu_;
        }

        //- Return the stoichiometric air/fuel ratio
        scalar stoicRatio() const
        {
            return stoicRatio_;
        }

        //- Return the fuel-oxidant equivalence ratio
        scalar Phi(const scalarFieldListSlice& Yu) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const uMulticomponentMixture&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
