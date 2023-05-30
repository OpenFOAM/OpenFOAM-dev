/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::hePsiThermo

Description
    Thermo implementation based on compressibility

SourceFiles
    hePsiThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hePsiThermo_H
#define hePsiThermo_H

#include "psiThermo.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class hePsiThermo Declaration
\*---------------------------------------------------------------------------*/

template<class HeThermo>
class hePsiThermo
:
    public HeThermo
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

public:

    // Constructors

        //- Construct from mesh and phase name
        hePsiThermo(const fvMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        hePsiThermo(const hePsiThermo<HeThermo>&) = delete;


    //- Destructor
    virtual ~hePsiThermo();


    // Member Functions

        //- Update properties
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const hePsiThermo<HeThermo>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "hePsiThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
