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
    Foam::CloudStateFieldRef

Description
    Wrapper for a reference to a dynamic Lagrangian field. It provides
    convenient access consistent with that of the derived fields so that
    expressions involving both kinds of field are readable. It also ensures
    that old-time fields are maintained for a sub-mesh for the duration of the
    inner correction iterations.

See also
    Foam::CloudStateField
    Foam::CloudDerivedField

SourceFiles
    CloudStateFieldRef.C

\*---------------------------------------------------------------------------*/

#ifndef CloudStateFieldRef_H
#define CloudStateFieldRef_H

#include "LagrangianModelRef.H"
#include "LagrangianFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class CloudStateFieldRef Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CloudStateFieldRef
{
private:

    // Private Data

        //- The field reference
        LagrangianDynamicField<Type>& ref_;

        //- The cached Lagrangian sub-all-field
        mutable autoPtr<LagrangianSubSubField<Type>> psiAllPtr_;

        //- The cached Lagrangian sub-field
        mutable autoPtr<LagrangianSubSubField<Type>> psiSubSubPtr_;

        //- The index of the sub-mesh
        mutable label psiSubSubMeshIndex_;


public:

    // Constructors

        //- Construct from a field reference
        CloudStateFieldRef(LagrangianDynamicField<Type>&);


    // Member Functions

        //- Access the field's name
        const word& name() const;

        //- Access the field's mesh
        const LagrangianMesh& mesh() const;

        //- Access the field's dimensions
        const dimensionSet& dimensions() const;

        //- Access a part of the field
        const LagrangianSubSubField<Type>& ref(const LagrangianSubMesh&) const;

        //- Access a part of the field
        LagrangianSubSubField<Type>& ref(const LagrangianSubMesh&);


    // Member Operators

        //- Cast to the field
        operator const LagrangianDynamicField<Type>&() const;

        //- Cast to the field
        operator LagrangianDynamicField<Type>&();

        //- Access the entire field
        const LagrangianDynamicField<Type>& operator()
        (
            const LagrangianMesh&
        ) const;

        //- Access a part of the field
        const LagrangianSubSubField<Type>& operator()
        (
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
    #include "CloudStateFieldRef.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
