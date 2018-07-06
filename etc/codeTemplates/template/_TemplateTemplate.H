/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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
    Foam::CLASSNAME

Description

SourceFiles
    CLASSNAMEI.H
    CLASSNAME.C
    CLASSNAMEIO.C

\*---------------------------------------------------------------------------*/

#ifndef CLASSNAME_H
#define CLASSNAME_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class someClass;


// Forward declaration of friend functions and operators

template<TemplateClassArgument>
class CLASSNAME;

template<TemplateClassArgument>
Istream& operator>>(Istream&, CLASSNAME<TemplateArgument>&);

template<TemplateClassArgument>
Ostream& operator<<(Ostream&, const CLASSNAME<TemplateArgument>&);


/*---------------------------------------------------------------------------*\
                         Class CLASSNAME Declaration
\*---------------------------------------------------------------------------*/

template<TemplateClassArgument>
class CLASSNAME
:
    public baseClassName
{
    // Private data

        dataType data_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        CLASSNAME(const CLASSNAME<TemplateArgument>&);

        //- Disallow default bitwise assignment
        void operator=(const CLASSNAME<TemplateArgument>&);


public:

    // Static data members

        //- Static data someStaticData
        static const dataType staticData;


    // Constructors

        //- Construct null
        CLASSNAME();

        //- Construct from components
        CLASSNAME(const dataType& data);

        //- Construct from Istream
        CLASSNAME(Istream&);

        //- Construct as copy
        CLASSNAME(const CLASSNAME<TemplateArgument>&);


    // Selectors

        //- Select null constructed
        static autoPtr<CLASSNAME<TemplateArgument>> New();


    //- Destructor
    ~CLASSNAME();


    // Member Functions

    // Member Operators

        void operator=(const CLASSNAME<TemplateArgument>&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>> <TemplateArgument>
        (Istream&, CLASSNAME<TemplateArgument>&);

        friend Ostream& operator<< <TemplateArgument>
        (Ostream&, const CLASSNAME<TemplateArgument>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CLASSNAMEI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CLASSNAME.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
