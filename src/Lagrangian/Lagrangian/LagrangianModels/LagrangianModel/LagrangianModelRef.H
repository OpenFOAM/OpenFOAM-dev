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
    Foam::LagrangianModelRef

Description
    Simple wrapper to provide an optional reference to a Lagrangian model

SourceFiles
    LagrangianModelRef.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianModelRef_H
#define LagrangianModelRef_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class LagrangianModel;

/*---------------------------------------------------------------------------*\
                     Class LagrangianModelRef Declaration
\*---------------------------------------------------------------------------*/

class LagrangianModelRef
{
    // Private Data

        //- Reference to the model
        const LagrangianModel& ref_;


public:

    // Constructors

        //- Construct null
        inline LagrangianModelRef();

        //- Construct for given reference
        inline LagrangianModelRef(const LagrangianModel& ref);


    // Member Functions

        //- Check if the reference is valid
        inline bool valid() const;


    // Member Operators

        //- Access the reference
        inline const LagrangianModel& operator()() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LagrangianModelRefI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
