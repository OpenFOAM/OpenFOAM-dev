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

\*---------------------------------------------------------------------------*/

#ifndef baseClass_H
#define baseClass_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class baseClass Declaration
\*---------------------------------------------------------------------------*/

template<class T1, class T2>
class baseClass
{
public:
    
    //- Runtime type information
    TypeName("baseClass");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        baseClass,
        dictionary,
        (
            word i, float f
        ),
        (i, f)
    );

    // Constructors

        // Construct from components
        baseClass (word i, float f) {};

        // Construct from copy
        baseClass (const baseClass&) = delete;

    // Selectors

        //- Return a reference to the selected baseClass
        static autoPtr<baseClass> New
        (
            word i, float f
        );


    //- Destructor
    virtual ~baseClass();

    // Member Functions

        //- Corrector
        virtual void correct() = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        baseClass& operator=(const baseClass&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "baseClass.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
