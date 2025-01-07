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
    Foam::LagrangianSource

Description
    Base class for Lagrangian sources. Minimal wrapper over LagrangianModel
    that provides an interface to the source value, and serves as a type that
    can be identified and cast to.

SourceFiles
    LagrangianSource.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianSource_H
#define LagrangianSource_H

#include "LagrangianModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class LagrangianSource Declaration
\*---------------------------------------------------------------------------*/

class LagrangianSource
:
    public LagrangianModel
{
public:

    //- Runtime type information
    TypeName("LagrangianSource");


    // Constructors

        //- Construct from components
        LagrangianSource(const word& name, const LagrangianMesh& mesh);


    //- Destructor
    virtual ~LagrangianSource();


    // Member Functions

        // Checks

            //- Return true if the LagrangianModel adds a source term to the
            //  given field's transport equation. Must be provided by sources.
            virtual bool addsSupToField(const word& fieldName) const = 0;


        // Sources

            //- Return the source value
            virtual tmp<LagrangianSubScalarField> source
            (
                const word& fieldName,
                const LagrangianSubMesh& subMesh
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
