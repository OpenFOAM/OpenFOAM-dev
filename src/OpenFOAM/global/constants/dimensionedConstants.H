/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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


class registerDimensionedConstant
:
    public simpleRegIOobject,
    public dimensionedScalar
{
    const word& group_;

    //- Lookup the dimensionedConstant from the registry
    void lookup();


public:

    registerDimensionedConstant
    (
        const word& group,
        const char* name,
        const dimensionSet& dimensions
    )
    :
        simpleRegIOobject(debug::addDimensionedConstantObject, name),
        dimensionedScalar(name, dimensions, 0),
        group_(group)
    {
        lookup();
    }

    virtual ~registerDimensionedConstant()
    {}

    virtual void readData(Istream& is)
    {
        lookup();
    }

    virtual void writeData(Ostream& os) const
    {
        os << static_cast<const dimensionedScalar&>(*this);
    }
};


//- Defined and register a dimensioned constant
#define defineDimensionedConstant(Group, Name, Dimensions)                     \
    registerDimensionedConstant register##Name##DimensionedConstant_           \
    (                                                                          \
        Group::group,                                                          \
        #Name,                                                                 \
        Dimensions                                                             \
    );                                                                         \
    const Foam::dimensionedScalar& Group::Name                                 \
    (                                                                          \
        register##Name##DimensionedConstant_                                   \
    );


class registerDimensionedConstantWithDefault
:
    public simpleRegIOobject,
    public dimensionedScalar
{
    const word& group_;
    const dimensionedScalar default_;

    //- Lookup the dimensionedConstant from the registry
    void lookup();


public:

    registerDimensionedConstantWithDefault
    (
        const word& group,
        const char* name,
        const dimensionedScalar& value
    )
    :
        simpleRegIOobject(debug::addDimensionedConstantObject, name),
        dimensionedScalar(name, value),
        group_(group),
        default_(name, value)
    {
        lookup();
    }

    virtual ~registerDimensionedConstantWithDefault()
    {}

    virtual void readData(Istream& is)
    {
        lookup();
    }

    virtual void writeData(Ostream& os) const
    {
        os << static_cast<const dimensionedScalar&>(*this);
    }
};


//- Defined and register a dimensioned constant with default
#define defineDimensionedConstantWithDefault(Group, Name, DefaultExpr)         \
    registerDimensionedConstantWithDefault                                     \
        register##Name##DimensionedConstant_                                   \
        (                                                                      \
            Group::group,                                                      \
            #Name,                                                             \
            DefaultExpr                                                        \
        );                                                                     \
                                                                               \
    const Foam::dimensionedScalar& Group::Name                                 \
    (                                                                          \
        register##Name##DimensionedConstant_                                   \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
