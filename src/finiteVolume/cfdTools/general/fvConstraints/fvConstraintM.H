/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#define DEFINE_FV_CONSTRAINT_CONSTRAIN(Type, nullArg)                          \
    virtual void constrain                                                     \
    (                                                                          \
        fvMatrix<Type>& eqn,                                                   \
        const word& fieldName                                                  \
    ) const;

#define IMPLEMENT_FV_CONSTRAINT_CONSTRAIN(Type, constraintType)                \
    void Foam::constraintType::constrain                                       \
    (                                                                          \
        fvMatrix<Type>& eqn,                                                   \
        const word& fieldName                                                  \
    ) const                                                                    \
    {                                                                          \
        constrainType(eqn, fieldName);                                         \
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define DEFINE_FV_CONSTRAINT_CONSTRAIN_FIELD(Type, nullArg)                    \
    virtual void constrain(VolField<Type>& field) const;

#define IMPLEMENT_FV_CONSTRAINT_CONSTRAIN_FIELD(Type, constraintType)          \
    void Foam::constraintType::constrain(VolField<Type>& field) const          \
    {                                                                          \
        constrainType(field);                                                  \
    }


// ************************************************************************* //
