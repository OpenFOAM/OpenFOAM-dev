/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::fvMeshToFvMesh

Description

SourceFiles
    fvMeshToFvMeshTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshToFvMesh_H
#define fvMeshToFvMesh_H

#include "meshToMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class fvMeshToFvMesh Declaration
\*---------------------------------------------------------------------------*/

class fvMeshToFvMesh
:
    public meshToMesh
{
public:

    //- Run-time type information
    TypeName("fvMeshToFvMesh");


    // Constructors

        //- Inherit base class' constructors
        using meshToMesh::meshToMesh;


    //- Destructor
    virtual ~fvMeshToFvMesh();


    // Member Functions

        // Evaluation

            // Source-to-target volume field mapping

                //- Interpolate a field with a defined operation. Values
                //  passed in via 'result' are used to initialise the return
                //  value.
                template<class Type>
                void mapSrcToTgt
                (
                    const VolField<Type>& field,
                    VolField<Type>& result
                ) const;

                //- Interpolate a field with a defined operation. The initial
                //  values of the result are set to zero.
                template<class Type>
                tmp<VolField<Type>> mapSrcToTgt
                (
                    const VolField<Type>& field
                ) const;


            // Source-to-target volume internal field mapping

                //- Interpolate a field with a defined operation. Values
                //  passed in via 'result' are used to initialise the return
                //  value.
                template<class Type>
                void mapSrcToTgt
                (
                    const VolInternalField<Type>& field,
                    VolInternalField<Type>& result
                ) const;

                //- Interpolate a field with a defined operation. The initial
                //  values of the result are set to zero.
                template<class Type>
                tmp<VolInternalField<Type>> mapSrcToTgt
                (
                    const VolInternalField<Type>& field
                ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshToFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
