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
    Foam::LagrangianModelRef

Description
    Simple wrapper to provide an optional reference to a Lagrangian model

SourceFiles
    LagrangianModelRef.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianModelRef_H
#define LagrangianModelRef_H

#include "LagrangianSubFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class LagrangianSource;
class LagrangianInjection;

/*---------------------------------------------------------------------------*\
                     Class LagrangianModelRef Declaration
\*---------------------------------------------------------------------------*/

class LagrangianModelRef
{
    // Private Data

        //- Reference to the source model
        const LagrangianSource& source_;

        //- Reference to the source field
        const LagrangianSubScalarField& S_;

        //- Reference to the injection model
        const LagrangianInjection& injection_;


public:

    // Constructors

        //- Construct null
        inline LagrangianModelRef();

        //- Construct for a source
        inline LagrangianModelRef
        (
            const LagrangianSource& source,
            const LagrangianSubScalarField& S
        );

        //- Construct for an injection
        inline LagrangianModelRef
        (
            const LagrangianInjection& injection
        );


    // Member Functions

        //- Check if the source reference is valid
        inline bool isSource() const;

        //- Access the source model
        inline const LagrangianSource& source() const;

        //- Access the source field
        inline const LagrangianSubScalarField& S() const;

        //- Check if the injection reference is valid
        inline bool isInjection() const;

        //- Access the injection model
        inline const LagrangianInjection& injection() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LagrangianModelRefI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
