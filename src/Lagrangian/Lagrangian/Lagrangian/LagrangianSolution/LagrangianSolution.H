/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::LagrangianSolution

Description
    Selector class for Lagrangian schemes

SourceFiles
    LagrangianSolution.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianSolution_H
#define LagrangianSolution_H

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class LagrangianSolution Declaration
\*---------------------------------------------------------------------------*/

class LagrangianSolution
:
    public IOdictionary
{
public:

    //- Type information
    ClassName("LagrangianSolution");


    // Constructors

        //- Construct for objectRegistry
        LagrangianSolution(const objectRegistry& db);

        //- Disallow default bitwise copy construction
        LagrangianSolution(const LagrangianSolution&) = delete;


    //- Destructor
    ~LagrangianSolution();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LagrangianSolution&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
