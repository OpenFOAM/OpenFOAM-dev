/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Description
    High performance macro functions for Field\<Type\> algebra.  These expand
    using either array element access (for vector machines) or pointer
    dereferencing for scalar machines as appropriate.

\*---------------------------------------------------------------------------*/

#ifndef FieldM_H
#define FieldM_H

#include "error.H"
#include "ListLoopM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef FULLDEBUG

template<class Type1, class Type2>
void checkFields
(
    const UList<Type1>& f1,
    const UList<Type2>& f2,
    const char* op
)
{
    if (f1.size() != f2.size())
    {
        FatalErrorInFunction
            << " Field<"<<pTraits<Type1>::typeName<<"> f1("<<f1.size()<<')'
            << " and Field<"<<pTraits<Type2>::typeName<<"> f2("<<f2.size()<<')'
            << endl << " for operation " << op
            << abort(FatalError);
    }
}

template<class Type1, class Type2, class Type3>
void checkFields
(
    const UList<Type1>& f1,
    const UList<Type2>& f2,
    const UList<Type3>& f3,
    const char* op
)
{
    if (f1.size() != f2.size() || f1.size() != f3.size())
    {
        FatalErrorInFunction
            << " Field<"<<pTraits<Type1>::typeName<<"> f1("<<f1.size()<<')'
            << ", Field<"<<pTraits<Type2>::typeName<<"> f2("<<f2.size()<<')'
            << " and Field<"<<pTraits<Type3>::typeName<<"> f3("<<f3.size()<<')'
            << endl << "    for operation " << op
            << abort(FatalError);
    }
}

#else

template<class Type1, class Type2>
void checkFields
(
    const UList<Type1>&,
    const UList<Type2>&,
    const char*
)
{}

template<class Type1, class Type2, class Type3>
void checkFields
(
    const UList<Type1>&,
    const UList<Type2>&,
    const UList<Type3>&,
    const char*
)
{}

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// member function : this f1 OP fUNC f2

#define TFOR_ALL_F_OP_FUNC_F(typeF1, f1, OP, FUNC, typeF2, f2)                 \
                                                                               \
    /* check the two fields have same Field<Type> mesh */                      \
    checkFields(f1, f2, "f1 " #OP " " #FUNC "(f2)");                           \
                                                                               \
    /* set access to f1, f2 and f3 at end of each field */                     \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* loop through fields performing f1 OP1 f2 OP2 f3 */                      \
    List_FOR_ALL(f1, i)                                                        \
        List_ELEM(f1, f1P, i) OP FUNC(List_ELEM(f2, f2P, i));                  \
    List_END_FOR_ALL                                                           \


#define TFOR_ALL_F_OP_F_FUNC(typeF1, f1, OP, typeF2, f2, FUNC)                 \
                                                                               \
    /* check the two fields have same Field<Type> mesh */                      \
    checkFields(f1, f2, "f1 " #OP " f2" #FUNC);                                \
                                                                               \
    /* set access to f1, f2 and f3 at end of each field */                     \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* loop through fields performing f1 OP1 f2 OP2 f3 */                      \
    List_FOR_ALL(f1, i)                                                        \
        List_ELEM(f1, f1P, i) OP List_ELEM(f2, f2P, i).FUNC();                 \
    List_END_FOR_ALL                                                           \


// member function : this field f1 OP fUNC f2, f3

#define TFOR_ALL_F_OP_FUNC_F_F(typeF1, f1, OP, FUNC, typeF2, f2, typeF3, f3)\
                                                                               \
    /* check the three fields have same Field<Type> mesh */                    \
    checkFields(f1, f2, f3, "f1 " #OP " " #FUNC "(f2, f3)");                   \
                                                                               \
    /* set access to f1, f2 and f3 at end of each field */                     \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
    List_CONST_ACCESS(typeF3, f3, f3P);                                        \
                                                                               \
    /* loop through fields performing f1 OP1 f2 OP2 f3 */                      \
    List_FOR_ALL(f1, i)                                                        \
        List_ELEM(f1, f1P, i)                                                  \
        OP FUNC(List_ELEM(f2, f2P, i), List_ELEM(f3, f3P, i));                 \
    List_END_FOR_ALL                                                           \


// member function : this field f1 OP fUNC f2, f3

#define TFOR_ALL_S_OP_FUNC_F_F(typeS, s, OP, FUNC, typeF1, f1, typeF2, f2)     \
                                                                               \
    /* check the two fields have same Field<Type> mesh */                      \
    checkFields(f1, f2, "s " #OP " " #FUNC "(f1, f2)");                        \
                                                                               \
    /* set access to f1 and f2 at end of each field */                         \
    List_CONST_ACCESS(typeF1, f1, f1P);                                        \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* loop through fields performing s OP FUNC(f1, f2) */                     \
    List_FOR_ALL(f1, i)                                                        \
        (s) OP FUNC(List_ELEM(f1, f1P, i), List_ELEM(f2, f2P, i));             \
    List_END_FOR_ALL                                                           \


// member function : this f1 OP fUNC f2, s

#define TFOR_ALL_F_OP_FUNC_F_S(typeF1, f1, OP, FUNC, typeF2, f2, typeS, s)     \
                                                                               \
    /* check the two fields have same Field<Type> mesh */                      \
    checkFields(f1, f2, "f1 " #OP " " #FUNC "(f2, s)");                        \
                                                                               \
    /* set access to f1, f2 and f3 at end of each field */                     \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* loop through fields performing f1 OP1 f2 OP2 f3 */                      \
    List_FOR_ALL(f1, i)                                                        \
        List_ELEM(f1, f1P, i) OP FUNC(List_ELEM(f2, f2P, i), (s));             \
    List_END_FOR_ALL


// member function : s1 OP fUNC f, s2

#define TFOR_ALL_S_OP_FUNC_F_S(typeS1, s1, OP, FUNC, typeF, f, typeS2, s2)     \
                                                                               \
    /* set access to f at end of field */                                      \
    List_CONST_ACCESS(typeF, f, fP);                                           \
                                                                               \
    /* loop through fields performing f1 OP1 f2 OP2 f3 */                      \
    List_FOR_ALL(f, i)                                                         \
        (s1) OP FUNC(List_ELEM(f, fP, i), (s2));                               \
    List_END_FOR_ALL                                                           \


// member function : this f1 OP fUNC s, f2

#define TFOR_ALL_F_OP_FUNC_S_F(typeF1, f1, OP, FUNC, typeS, s, typeF2, f2)     \
                                                                               \
    /* check the two fields have same Field<Type> mesh */                      \
    checkFields(f1, f2, "f1 " #OP " " #FUNC "(s, f2)");                        \
                                                                               \
    /* set access to f1, f2 and f3 at end of each field */                     \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* loop through fields performing f1 OP1 f2 OP2 f3 */                      \
    List_FOR_ALL(f1, i)                                                        \
        List_ELEM(f1, f1P, i) OP FUNC((s), List_ELEM(f2, f2P, i));             \
    List_END_FOR_ALL                                                           \


// member function : this f1 OP fUNC s, f2

#define TFOR_ALL_F_OP_FUNC_S_S(typeF1, f1, OP, FUNC, typeS1, s1, typeS2, s2)\
                                                                               \
    /* set access to f1 at end of field */                                     \
    List_ACCESS(typeF1, f1, f1P);                                              \
                                                                               \
    /* loop through fields performing f1 OP1 FUNC(s1, s2) */                   \
    List_FOR_ALL(f1, i)                                                        \
        List_ELEM(f1, f1P, i) OP FUNC((s1), (s2));                             \
    List_END_FOR_ALL                                                           \


// member function : this f1 OP1 f2 OP2 FUNC s

#define TFOR_ALL_F_OP_F_FUNC_S(typeF1, f1, OP, typeF2, f2, FUNC, typeS, s)     \
                                                                               \
    /* check the two fields have same Field<Type> mesh */                      \
    checkFields(f1, f2, "f1 " #OP " f2 " #FUNC "(s)");                         \
                                                                               \
    /* set access to f1, f2 and f3 at end of each field */                     \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* loop through fields performing f1 OP1 f2 OP2 f3 */                      \
    List_FOR_ALL(f1, i)                                                        \
        List_ELEM(f1, f1P, i) OP List_ELEM(f2, f2P, i) FUNC((s));              \
    List_END_FOR_ALL                                                           \


// define high performance macro functions for Field<Type> operations

// member operator : this field f1 OP1 f2 OP2 f3

#define TFOR_ALL_F_OP_F_OP_F(typeF1, f1, OP1, typeF2, f2, OP2, typeF3, f3)     \
                                                                               \
    /* check the three fields have same Field<Type> mesh */                    \
    checkFields(f1, f2, f3, "f1 " #OP1 " f2 " #OP2 " f3");                     \
                                                                               \
    /* set access to f1, f2 and f3 at end of each field */                     \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
    List_CONST_ACCESS(typeF3, f3, f3P);                                        \
                                                                               \
    /* loop through fields performing f1 OP1 f2 OP2 f3 */                      \
    List_FOR_ALL(f1, i)                                                        \
        List_ELEM(f1, f1P, i) OP1 List_ELEM(f2, f2P, i)                        \
                              OP2 List_ELEM(f3, f3P, i);                       \
    List_END_FOR_ALL                                                           \


// member operator : this field f1 OP1 s OP2 f2

#define TFOR_ALL_F_OP_S_OP_F(typeF1, f1, OP1, typeS, s, OP2, typeF2, f2)       \
                                                                               \
    /* check the two fields have same Field<Type> mesh */                      \
    checkFields(f1, f2, "f1 " #OP1 " s " #OP2 " f2");                          \
                                                                               \
    /* set access to f1 and f2 at end of each field */                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* loop through fields performing f1 OP1 s OP2 f2 */                       \
    List_FOR_ALL(f1, i)                                                        \
        List_ELEM(f1, f1P, i) OP1 (s) OP2 List_ELEM(f2, f2P, i);               \
    List_END_FOR_ALL                                                           \


// member operator : this field f1 OP1 f2 OP2 s

#define TFOR_ALL_F_OP_F_OP_S(typeF1, f1, OP1, typeF2, f2, OP2, typeS, s)       \
                                                                               \
    /* check the two fields have same Field<Type> mesh */                      \
    checkFields(f1, f2, "f1 " #OP1 " f2 " #OP2 " s");                          \
                                                                               \
    /* set access to f1 and f2 at end of each field */                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* loop through fields performing f1 OP1 s OP2 f2 */                       \
    List_FOR_ALL(f1, i)                                                        \
        List_ELEM(f1, f1P, i) OP1 List_ELEM(f2, f2P, i) OP2 (s);               \
    List_END_FOR_ALL                                                           \


// member operator : this field f1 OP f2

#define TFOR_ALL_F_OP_F(typeF1, f1, OP, typeF2, f2)                            \
                                                                               \
    /* check the two fields have same Field<Type> mesh */                      \
    checkFields(f1, f2, "f1 " #OP " f2");                                      \
                                                                               \
    /* set pointer to f1P at end of f1 and */                                  \
    /* f2.p at end of f2 */                                                    \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* loop through fields performing f1 OP f2 */                              \
    List_FOR_ALL(f1, i)                                                        \
        List_ELEM(f1, f1P, i) OP List_ELEM(f2, f2P, i);                        \
    List_END_FOR_ALL                                                           \

// member operator : this field f1 OP1 OP2 f2

#define TFOR_ALL_F_OP_OP_F(typeF1, f1, OP1, OP2, typeF2, f2)                   \
                                                                               \
    /* check the two fields have same Field<Type> mesh */                      \
    checkFields(f1, f2, #OP1 " " #OP2 " f2");                                  \
                                                                               \
    /* set pointer to f1P at end of f1 and */                                  \
    /* f2.p at end of f2 */                                                    \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* loop through fields performing f1 OP1 OP2 f2 */                         \
    List_FOR_ALL(f1, i)                                                        \
        List_ELEM(f1, f1P, i) OP1 OP2 List_ELEM(f2, f2P, i);                   \
    List_END_FOR_ALL                                                           \


// member operator : this field f OP s

#define TFOR_ALL_F_OP_S(typeF, f, OP, typeS, s)                                \
                                                                               \
    /* set access to f at end of field */                                      \
    List_ACCESS(typeF, f, fP);                                                 \
                                                                               \
    /* loop through field performing f OP s */                                 \
    List_FOR_ALL(f, i)                                                         \
        List_ELEM(f, fP, i) OP (s);                                            \
    List_END_FOR_ALL                                                           \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// define high performance macro functions for Field<Type> friend functions

// friend operator function : s OP f, allocates storage for s

#define TFOR_ALL_S_OP_F(typeS, s, OP, typeF, f)                                \
                                                                               \
    /* set access to f at end of field */                                      \
    List_CONST_ACCESS(typeF, f, fP);                                           \
                                                                               \
    /* loop through field performing s OP f */                                 \
    List_FOR_ALL(f, i)                                                         \
        (s) OP List_ELEM(f, fP, i);                                            \
    List_END_FOR_ALL


// friend operator function : s OP1 f1 OP2 f2, allocates storage for s

#define TFOR_ALL_S_OP_F_OP_F(typeS, s, OP1, typeF1, f1, OP2, typeF2, f2)       \
                                                                               \
    /* set access to f1 and f2 at end of each field */                         \
    List_CONST_ACCESS(typeF1, f1, f1P);                                        \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* loop through field performing s OP f */                                 \
    List_FOR_ALL(f1, i)                                                        \
        (s) OP1 List_ELEM(f1, f1P, i) OP2 List_ELEM(f2, f2P, i);               \
    List_END_FOR_ALL


// friend operator function : s OP FUNC(f), allocates storage for s

#define TFOR_ALL_S_OP_FUNC_F(typeS, s, OP, FUNC, typeF, f)                     \
                                                                               \
    /* set access to f at end of field */                                      \
    List_CONST_ACCESS(typeF, f, fP);                                           \
                                                                               \
    /* loop through field performing s OP f */                                 \
    List_FOR_ALL(f, i)                                                         \
        (s) OP FUNC(List_ELEM(f, fP, i));                                      \
    List_END_FOR_ALL


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
