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

Class
    Foam::UniformDimensionedField

Description
    Dimensioned<Type> registered with the database as a registered IOobject
    which has the functionality of a uniform field and allows values from the
    top-level code to be passed to boundary conditions etc.

SourceFiles
    UniformDimensionedField.C

\*---------------------------------------------------------------------------*/

#ifndef UniformDimensionedField_H
#define UniformDimensionedField_H

#include "regIOobject.H"
#include "dimensionedType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class UniformDimensionedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class UniformDimensionedField
:
    public regIOobject,
    public dimensioned<Type>
{

public:

    //- Runtime type information
    TypeName("UniformDimensionedField");


    // Constructors

        //- Construct from components. Either reads or uses supplied value.
        UniformDimensionedField(const IOobject&, const dimensioned<Type>&);

        //- Construct as copy
        UniformDimensionedField(const UniformDimensionedField<Type>&);

        //- Construct from Istream
        UniformDimensionedField(const IOobject&);


    //- Destructor
    virtual ~UniformDimensionedField();


    // Member Functions

        //- Name function provided to resolve the ambiguity between the
        //  name functions in regIOobject and dimensioned<Type>
        virtual const word& name() const
        {
            return dimensioned<Type>::name();
        }

        bool writeData(Ostream&) const;

        //- Is object global
        virtual bool global() const
        {
            return true;
        }

        //- Return complete path + object name if the file exists
        //  either in the case/processor or case otherwise null
        virtual fileName filePath() const
        {
            return globalFilePath(type());
        }


    // Member Operators

        void operator=(const UniformDimensionedField<Type>&);
        void operator=(const dimensioned<Type>&);

        Type operator[](const label) const
        {
            return this->value();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UniformDimensionedField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
