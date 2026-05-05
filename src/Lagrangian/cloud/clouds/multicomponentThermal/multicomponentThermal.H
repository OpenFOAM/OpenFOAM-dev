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
    Foam::clouds::multicomponentThermal

Description
    Base class for clouds with multicomponent thermodynamic modelling

SourceFiles
    multicomponentThermal.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentThermal_H
#define multicomponentThermal_H

#include "thermal.H"
#include "multicomponentLagrangianThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{

/*---------------------------------------------------------------------------*\
                     Class multicomponentThermal Declaration
\*---------------------------------------------------------------------------*/

class multicomponentThermal
:
    public Thermal<multicomponentLagrangianThermo>
{
public:

    // Public Static Data

        //- Run-time type information
        TypeName("multicomponentThermal");


    // Constructors

        //- Construct from a reference to the cloud and its base classes
        multicomponentThermal(const cloud&, const shaped&, const carried&);


    // Public Data

        //- Mass fractions
        PtrList<CloudStateFieldRef<scalar>> Y;


    //- Destructor
    virtual ~multicomponentThermal();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
