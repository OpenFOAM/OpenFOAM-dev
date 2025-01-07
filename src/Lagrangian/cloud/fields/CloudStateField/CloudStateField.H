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
    Foam::CloudStateField

Description
    A field which is stored as part of the state of the cloud. This is a light
    wrapper around a dynamic Lagrangian field. It provides convenient access
    consistent with that of the derived fields so that expressions involving
    both kinds of field are readable.

See also
    Foam::CloudDerivedField

SourceFiles
    CloudStateField.C

\*---------------------------------------------------------------------------*/

#ifndef CloudStateField_H
#define CloudStateField_H

#include "LagrangianModelRef.H"
#include "LagrangianFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class CloudStateField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CloudStateField
:
    public LagrangianDynamicField<Type>
{
private:

    // Private Data

        //- The cached Lagrangian sub-field
        mutable autoPtr<LagrangianSubSubField<Type>> psiSubSubPtr_;

        //- The index of the sub-mesh
        mutable label psiSubSubMeshIndex_;


public:

    // Constructors

        //- Use the base field's constructors
        using LagrangianDynamicField<Type>::LagrangianDynamicField;


    // Member Functions

        //- Access a part of the field
        LagrangianSubSubField<Type>& ref(const LagrangianSubMesh& subMesh);

        //- Clear
        void clear();


    // Member Operators

        //- Access the entire field
        const LagrangianDynamicField<Type>& operator()
        (
            const LagrangianMesh&
        ) const;

        //- Access a part of the field
        tmp<LagrangianSubSubField<Type>> operator()
        (
            const LagrangianSubMesh& subMesh
        ) const;

        //- Compute and return part of the field associated with a model
        tmp<LagrangianSubField<Type>> operator()
        (
            const LagrangianModel& model,
            const LagrangianSubMesh& subMesh
        ) const;

        //- Compute and return part of the field associated with a model
        tmp<LagrangianSubSubField<Type>> operator()
        (
            const LagrangianModelRef& model,
            const LagrangianSubMesh& subMesh
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CloudStateField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
