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

\*---------------------------------------------------------------------------*/

#include "expressions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expression
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Inherit expressionability from the underlying Field
template<class Type>
struct Expressionable<fvsPatchField<Type>>
:
    public Expressionable<Field<Type>>
{};

//- Inherit expressionability from the base Field
template<class Type>
struct Expressionable<fvsPatchField<Type>, Arguments<label>>
:
    public Expressionable<Field<Type>, Arguments<label>>
{};

//- Inherit expressionability from the base Field
template<class Type>
struct Expressionable
<
    fvsPatchField<Type>,
    Arguments<typename Field<Type>::const_iterator>
>
:
    public Expressionable
    <
        Field<Type>,
        Arguments<typename Field<Type>::const_iterator>
    >
{};

//- Inherit access by index from the base Field
template<class Type>
struct Operate<fvsPatchField<Type>, Arguments<label>>
:
    public Operate<Field<Type>, Arguments<label>>
{};

//- Inherit access to the start and end iterators from the base Field
template<class Type, bool Index>
struct Operate<fvsPatchField<Type>, Arguments<Bound<Index>>>
:
    public Operate<Field<Type>, Arguments<Bound<Index>>>
{};

//- Inherit access to a value from the base Field
template<class Type, class ... Args>
struct Operate
<
    fvsPatchField<Type>,
    Arguments<typename Field<Type>::const_iterator, Args ...>
>
:
    public Operate
    <
        Field<Type>,
        Arguments<typename Field<Type>::const_iterator, Args ...>
    >
{};

//- Define access to the base expression of a fvsPatchField as the base Field
template<class Type>
struct Operate<fvsPatchField<Type>, Arguments<Base>>
{
    using type = const Field<Type>&;

    static type access(const fvsPatchField<Type>& f, const Base&)
    {
        return static_cast<type>(f);
    }
};

//- Define the types that result from operating on a fvsPatchField as the same
//  as that result from the base Field
template<class Type>
struct Result<fvsPatchField<Type>>
{
    template<class Expression>
    using operateType =
        typename Operate<std::decay_t<Expression>, Arguments<Base>>::type;

    template<class Expression>
    using elementType =
        typename Result
        <
            Field<Type>
        >::template elementType<operateType<Expression>>;

    template<class Expression>
    using type =
        typename Result
        <
            Field<Type>
        >::template type<operateType<Expression>>;
};

//- Define reduction of a fvsPatchField as the same as that of the base Field
template<class Type>
struct Reduce<fvsPatchField<Type>>
:
    public Reduce<Field<Type>>
{};

//- Implement selection of a fvsPatchField within an expression. Construct a
//  calculated field type on non-constrained patch types, and pass through a
//  reference to the internal field if this was provided.
template<class Type>
struct Select
<
    fvsPatchField<Type>,
    true,
    EnableIfExpressionable<fvsPatchField<Type>>
>
{
    using type = tmp<fvsPatchField<Type>>;

    template<class Expression>
    static type New
    (
        const Expression& e,
        const DimensionedField<Type, fvMesh, Field>& iF
    )
    {
        tmp<fvsPatchField<Type>> tResult =
            fvsPatchField<Type>::New
            (
                fvsPatchField<Type>::calculatedType(),
                getAll<Patch<fvPatch>>(e),
                iF
            );

        tResult.ref() = e;

        return tResult;
    }

    template<class Expression>
    static type New(const Expression& e)
    {
        return
            New
            (
                e,
                NullObjectRef<DimensionedField<Type, fvMesh, Field>>()
            );
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// ************************************************************************* //
