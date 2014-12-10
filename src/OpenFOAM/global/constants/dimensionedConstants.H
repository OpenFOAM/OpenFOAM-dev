/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Global
    dimensionedConstants

Description
    Dictionary reading and supplying the dimensioned constants used within
    OpenFOAM, particularly for thermodynamics.

    The values are read from the OpenFOAM etc/controlDict and should be
    changed to run with a different set of units from the default SI units.

SourceFiles
    dimensionedConstants.C

\*---------------------------------------------------------------------------*/

#ifndef dimensionedConstants_H
#define dimensionedConstants_H

#include "dictionary.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionary& dimensionedConstants();


dimensionedScalar dimensionedConstant(const word& group, const word& varName);


template<class T>
T dimensionedConstant
(
    const word& group,
    const word& varName,
    const T& defaultValue
)
{
    dictionary& dict = dimensionedConstants();

    const word unitSet(dict.lookup("unitSet"));

    dictionary& unitDict(dict.subDict(unitSet + "Coeffs"));

    if (unitDict.found(group))
    {
        dictionary& groupDict = unitDict.subDict(group);
        if (groupDict.found(varName))
        {
            return pTraits<T>(groupDict.lookup(varName));
        }
        else
        {
            groupDict.add(varName, defaultValue);
            return defaultValue;
        }
    }
    else
    {
        unitDict.add(group, dictionary::null);
        unitDict.subDict(group).add(varName, defaultValue);

        return defaultValue;
    }
}


//- Defined dimensioned constant , lookup as \a Name
#define defineDimensionedConstant(Group,Switch,Tag,Name)                      \
    const Foam::dimensionedScalar Switch;                                     \
    class add##Tag##ToDimensionedConstant                                     \
    :                                                                         \
        public Foam::simpleRegIOobject                                        \
    {                                                                         \
    public:                                                                   \
        add##Tag##ToDimensionedConstant(const char* name)                     \
        :                                                                     \
            Foam::simpleRegIOobject                                           \
            (Foam::debug::addDimensionedConstantObject,name)                  \
        {                                                                     \
            Foam::dimensionedScalar ds                                        \
            (                                                                 \
                Foam::dimensionedConstant                                     \
                (                                                             \
                    Group,                                                    \
                    Name                                                      \
                )                                                             \
            );                                                                \
            Foam::dimensionedScalar& s = const_cast<Foam::dimensionedScalar&> \
            (                                                                 \
                Switch                                                        \
            );                                                                \
            s.dimensions().reset(ds.dimensions());                            \
            s = ds;                                                           \
        }                                                                     \
        virtual ~add##Tag##ToDimensionedConstant()                            \
        {}                                                                    \
        virtual void readData(Foam::Istream& is)                              \
        {                                                                     \
            const_cast<Foam::dimensionedScalar&>(Switch) =                    \
                Foam::dimensionedConstant                                     \
                (                                                             \
                    Group,                                                    \
                    Name                                                      \
                );                                                            \
        }                                                                     \
        virtual void writeData(Foam::Ostream& os) const                       \
        {                                                                     \
            os << Switch;                                                     \
        }                                                                     \
    };                                                                        \
    add##Tag##ToDimensionedConstant add##Tag##ToDimensionedConstant_(Name)


//- Defined dimensioned constant with default , lookup as \a Name
#define defineDimensionedConstantWithDefault\
(Group,Switch,DefaultExpr,Tag,Name)                                           \
    const Foam::dimensionedScalar Switch;                                     \
    class add##Tag##ToDimensionedConstantWithDefault                          \
    :                                                                         \
        public Foam::simpleRegIOobject                                        \
    {                                                                         \
    public:                                                                   \
        add##Tag##ToDimensionedConstantWithDefault(const char* name)          \
        :                                                                     \
            Foam::simpleRegIOobject                                           \
            (Foam::debug::addDimensionedConstantObject,name)                  \
        {                                                                     \
            Foam::dimensionedScalar ds                                        \
            (                                                                 \
                Foam::dimensionedConstant                                     \
                (                                                             \
                    Group,                                                    \
                    Name,                                                     \
                    Foam::dimensionedScalar(Name,DefaultExpr)                 \
                )                                                             \
            );                                                                \
            Foam::dimensionedScalar& s = const_cast<Foam::dimensionedScalar&> \
            (                                                                 \
                Switch                                                        \
            );                                                                \
            s.dimensions().reset(ds.dimensions());                            \
            s = ds;                                                           \
        }                                                                     \
        virtual ~add##Tag##ToDimensionedConstantWithDefault()                 \
        {}                                                                    \
        virtual void readData(Foam::Istream& is)                              \
        {                                                                     \
            const_cast<Foam::dimensionedScalar&>(Switch) =                    \
                Foam::dimensionedConstant                                     \
                (                                                             \
                    Group,                                                    \
                    Name,                                                     \
                    Foam::dimensionedScalar(Name,DefaultExpr)                 \
                );                                                            \
        }                                                                     \
        virtual void writeData(Foam::Ostream& os) const                       \
        {                                                                     \
            os << Switch;                                                     \
        }                                                                     \
    };                                                                        \
    add##Tag##ToDimensionedConstantWithDefault                                \
        add##Tag##ToDimensionedConstantWithDefault_(Name)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
