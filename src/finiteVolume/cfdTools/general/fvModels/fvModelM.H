/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#define DEFINE_FV_MODEL_ADD_FIELD_SUP(Type, nullArg)                           \
    virtual void addSup                                                        \
    (                                                                          \
        const VolField<Type>& field,                                           \
        fvMatrix<Type>& eqn                                                    \
    ) const;

#define IMPLEMENT_FV_MODEL_ADD_FIELD_SUP(Type, modelType)                      \
    void Foam::modelType::addSup                                               \
    (                                                                          \
        const VolField<Type>& field,                                           \
        fvMatrix<Type>& eqn                                                    \
    ) const                                                                    \
    {                                                                          \
        addSupType(field, eqn);                                                \
    }

#define DEFINE_FV_MODEL_ADD_RHO_FIELD_SUP(Type, nullArg)                       \
    virtual void addSup                                                        \
    (                                                                          \
        const volScalarField& rho,                                             \
        const VolField<Type>& field,                                           \
        fvMatrix<Type>& eqn                                                    \
    ) const;

#define IMPLEMENT_FV_MODEL_ADD_RHO_FIELD_SUP(Type, modelType)                  \
    void Foam::modelType::addSup                                               \
    (                                                                          \
        const volScalarField& rho,                                             \
        const VolField<Type>& field,                                           \
        fvMatrix<Type>& eqn                                                    \
    ) const                                                                    \
    {                                                                          \
        addSupType(rho, field, eqn);                                           \
    }

#define DEFINE_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP(Type, nullArg)                 \
    virtual void addSup                                                        \
    (                                                                          \
        const volScalarField& alpha,                                           \
        const volScalarField& rho,                                             \
        const VolField<Type>& field,                                           \
        fvMatrix<Type>& eqn                                                    \
    ) const;

#define IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP(Type, modelType)            \
    void Foam::modelType::addSup                                               \
    (                                                                          \
        const volScalarField& alpha,                                           \
        const volScalarField& rho,                                             \
        const VolField<Type>& field,                                           \
        fvMatrix<Type>& eqn                                                    \
    ) const                                                                    \
    {                                                                          \
        addSupType(alpha, rho, field, eqn);                                    \
    }


// ************************************************************************* //
