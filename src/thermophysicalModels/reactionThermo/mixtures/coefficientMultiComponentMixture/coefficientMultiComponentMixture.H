/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::coefficientMultiComponentMixture

Description
    Foam::coefficientMultiComponentMixture

SourceFiles
    coefficientMultiComponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef coefficientMultiComponentMixture_H
#define coefficientMultiComponentMixture_H

#include "multiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class coefficientMultiComponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class coefficientMultiComponentMixture
:
    public multiComponentMixture<ThermoType>
{

public:

    //- The type of mixture
    typedef ThermoType mixtureType;


private:

    // Private Data

        //- Temporary storage for the cell/face mixture thermo data
        mutable mixtureType mixture_;


public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        coefficientMultiComponentMixture
        (
            const dictionary&,
            const fvMesh&,
            const word&
        );

        //- Disallow default bitwise copy construction
        coefficientMultiComponentMixture
        (
            const coefficientMultiComponentMixture<ThermoType>&
        ) = delete;


    //- Destructor
    virtual ~coefficientMultiComponentMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "multiComponentMixture<" + ThermoType::typeName() + '>';
        }

        const mixtureType& cellMixture(const label celli) const;

        const mixtureType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "coefficientMultiComponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
