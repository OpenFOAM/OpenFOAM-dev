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
    Foam::heLiquidThermo

Description
    Liquid thermo implementation

SourceFiles
    heLiquidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef heLiquidThermo_H
#define heLiquidThermo_H

#include "liquidThermo.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heLiquidThermo Declaration
\*---------------------------------------------------------------------------*/

template<class HeThermo>
class heLiquidThermo
:
    public HeThermo
{
public:

    // Constructors

        //- Construct from mesh and phase name
        heLiquidThermo(const fvMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        heLiquidThermo(const heLiquidThermo<HeThermo>&) = delete;


    //- Destructor
    virtual ~heLiquidThermo();


    // Member Functions

        //- Surface tension [N/m]
        virtual tmp<volScalarField> sigma() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const heLiquidThermo<HeThermo>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "heLiquidThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
