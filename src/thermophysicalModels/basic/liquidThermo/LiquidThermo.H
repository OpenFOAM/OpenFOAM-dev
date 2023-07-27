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
    Foam::LiquidThermo

Description
    Liquid thermo implementation

SourceFiles
    LiquidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef LiquidThermo_H
#define LiquidThermo_H

#include "BasicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class LiquidThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BaseThermo>
class LiquidThermo
:
    public BaseThermo
{
public:

    // Constructors

        //- Construct from mesh and phase name
        LiquidThermo(const fvMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        LiquidThermo(const LiquidThermo<BaseThermo>&) = delete;


    //- Destructor
    virtual ~LiquidThermo();


    // Member Functions

        //- Surface tension [N/m]
        virtual tmp<volScalarField> sigma() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LiquidThermo<BaseThermo>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LiquidThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
