/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#define DEFINE_FV_OPTION_ADD_SUP(Type, nullArg)                                \
    virtual void addSup                                                        \
    (                                                                          \
        fvMatrix<Type>& eqn,                                                   \
        const word& fieldName                                                  \
    ) const;

#define IMPLEMENT_FV_OPTION_ADD_SUP(Type, optionType)                          \
    void Foam::fv::optionType::addSup                                          \
    (                                                                          \
        fvMatrix<Type>& eqn,                                                   \
        const word& fieldName                                                  \
    ) const                                                                    \
    {                                                                          \
        addSupType(eqn, fieldName);                                            \
    }

#define DEFINE_FV_OPTION_ADD_RHO_SUP(Type, nullArg)                            \
    virtual void addSup                                                        \
    (                                                                          \
        const volScalarField& rho,                                             \
        fvMatrix<Type>& eqn,                                                   \
        const word& fieldName                                                  \
    ) const;

#define IMPLEMENT_FV_OPTION_ADD_RHO_SUP(Type, optionType)                      \
    void Foam::fv::optionType::addSup                                          \
    (                                                                          \
        const volScalarField& rho,                                             \
        fvMatrix<Type>& eqn,                                                   \
        const word& fieldName                                                  \
    ) const                                                                    \
    {                                                                          \
        addSupType(rho, eqn, fieldName);                                       \
    }

#define DEFINE_FV_OPTION_ADD_ALPHA_RHO_SUP(Type, nullArg)                      \
    virtual void addSup                                                        \
    (                                                                          \
        const volScalarField& alpha,                                           \
        const volScalarField& rho,                                             \
        fvMatrix<Type>& eqn,                                                   \
        const word& fieldName                                                  \
    ) const;

#define IMPLEMENT_FV_OPTION_ADD_ALPHA_RHO_SUP(Type, optionType)                \
    void Foam::fv::optionType::addSup                                          \
    (                                                                          \
        const volScalarField& alpha,                                           \
        const volScalarField& rho,                                             \
        fvMatrix<Type>& eqn,                                                   \
        const word& fieldName                                                  \
    ) const                                                                    \
    {                                                                          \
        addSupType(alpha, rho, eqn, fieldName);                                \
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define DEFINE_FV_OPTION_CONSTRAIN(Type, nullArg)                              \
    virtual void constrain                                                     \
    (                                                                          \
        fvMatrix<Type>& eqn,                                                   \
        const word& fieldName                                                  \
    ) const;

#define IMPLEMENT_FV_OPTION_CONSTRAIN(Type, optionType)                        \
    void Foam::fv::optionType::constrain                                       \
    (                                                                          \
        fvMatrix<Type>& eqn,                                                   \
        const word& fieldName                                                  \
    ) const                                                                    \
    {                                                                          \
        constrainType(eqn, fieldName);                                         \
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define DEFINE_FV_OPTION_CORRECT(Type, nullArg)                                \
    virtual void correct(VolField<Type>& field) const;

#define IMPLEMENT_FV_OPTION_CORRECT(Type, optionType)                          \
    void Foam::fv::optionType::correct(VolField<Type>& field) const            \
    {                                                                          \
        correctType(field);                                                    \
    }

// ************************************************************************* //
