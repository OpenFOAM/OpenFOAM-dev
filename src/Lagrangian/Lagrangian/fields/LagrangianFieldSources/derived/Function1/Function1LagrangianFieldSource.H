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
    Foam::Function1LagrangianFieldSource

Description
    Mix-in for source conditions that provides functions for evaluating
    Function1s at variable times.

SourceFiles
    Function1LagrangianFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef Function1LagrangianFieldSource_H
#define Function1LagrangianFieldSource_H

#include "LagrangianFieldsFwd.H"
#include "LagrangianSubFieldsFwd.H"
#include "dimensionSet.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class LagrangianModel;

/*---------------------------------------------------------------------------*\
               Class Function1LagrangianFieldSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Function1LagrangianFieldSource
{
    // Private Member Data

        //- Reference to the derived field source
        const LagrangianFieldSource<Type>& field_;


protected:

    // Protected Member Functions

        //- Return the source value
        template<class OtherType>
        static tmp<LagrangianSubField<OtherType>> value
        (
            const LagrangianModel& model,
            const LagrangianSubMesh& subMesh,
            const dimensionSet& dims,
            const Function1<OtherType>& function
        );

        //- Return the source value
        tmp<LagrangianSubField<Type>> value
        (
            const LagrangianModel& model,
            const LagrangianSubMesh& subMesh,
            const Function1<Type>& function
        ) const;


public:

    // Constructors

        //- Construct with reference to the derived field source
        template<class Derived>
        Function1LagrangianFieldSource(const Derived&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Function1LagrangianFieldSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
