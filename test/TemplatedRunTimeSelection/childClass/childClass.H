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

SourceFiles
    childClass.C

\*---------------------------------------------------------------------------*/

#ifndef childClass_H
#define childClass_H

#include "baseClass.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace models
{

/*---------------------------------------------------------------------------*\
                   Class childClass Declaration
\*---------------------------------------------------------------------------*/

template<class T1, class T2>
class childClass
:
    public baseClass<T1, T2>
{
public:

    //- Runtime type information
    ClassName("childClass");

    inline void spec();

    // Constructors

        // Construct from components
        childClass
        (
            word i,
            float f
        );

        // Construct from copy
        childClass (const childClass&) = delete;

	// Destructor
	virtual ~childClass();

    // Member Functions

		//- Phase fields update
		virtual void correct() override { spec(); }

    // Member operators

        //- Disallow default bitwise assignment
        childClass& operator=(const childClass&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phases

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "childClass.C"
    #include "childClassImpl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
