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
    Foam::SolidSurfaceThermo

Description
    Solid Lagrangian surface thermodynamic model implementation

SourceFiles
    SolidSurfaceThermo.C

\*---------------------------------------------------------------------------*/

#ifndef SolidSurfaceThermo_H
#define SolidSurfaceThermo_H

#include "basicLagrangianSubThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class SolidSurfaceThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BaseThermo>
class SolidSurfaceThermo
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
        SolidSurfaceThermo(const basicLagrangianSubThermo& farThermo);


    //- Destructor
    virtual ~SolidSurfaceThermo();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SolidSurfaceThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
