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
    Foam::FluidSurfaceThermo

Description
    Fluid Lagrangian surface thermodynamic model implementation

SourceFiles
    FluidSurfaceThermo.C

\*---------------------------------------------------------------------------*/

#ifndef FluidSurfaceThermo_H
#define FluidSurfaceThermo_H

#include "basicLagrangianSubThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class FluidSurfaceThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BaseThermo>
class FluidSurfaceThermo
:
    public BaseThermo
{
protected:

    // Protected Member Functions

        //- Correct surface properties
        void correct(const LagrangianSubMesh&) const;


public:

    // Constructors

        //- Construct from a reference to the far thermo
        FluidSurfaceThermo(const basicLagrangianSubThermo& farThermo);


    //- Destructor
    virtual ~FluidSurfaceThermo();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FluidSurfaceThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
