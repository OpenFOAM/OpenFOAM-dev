/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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
    both kinds of field are readable. It also ensures that old-time fields are
    maintained for a sub-mesh for the duration of the inner correction
    iterations.

See also
    Foam::CloudDerivedField

SourceFiles
    CloudStateField.C

\*---------------------------------------------------------------------------*/

#ifndef CloudStateField_H
#define CloudStateField_H

#include "CloudStateFieldRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class CloudStateField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CloudStateField
:
    public LagrangianDynamicField<Type>,
    public CloudStateFieldRef<Type>
{
public:

    // Constructors

        //- Forward construction arguments to the base class
        template<class ... Args>
        CloudStateField(const Args& ... args)
        :
            LagrangianDynamicField<Type>(args ...),
            CloudStateFieldRef<Type>
            (
                static_cast<LagrangianDynamicField<Type>&>(*this)
            )
        {}


    // Member Functions

        //- Access the name via the geometric field
        using LagrangianDynamicField<Type>::name;

        //- Access the name via the geometric field
        using LagrangianDynamicField<Type>::mesh;

        //- Access the name via the geometric field
        using LagrangianDynamicField<Type>::dimensions;


    // Member Operators

        //- Use the call operator(s) from the state field reference and not
        //  from the geometric field
        using CloudStateFieldRef<Type>::operator();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
