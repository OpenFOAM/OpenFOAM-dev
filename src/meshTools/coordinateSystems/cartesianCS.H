/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::cartesianCS

Description
    Cylindrical coordinate system

SourceFiles
    cartesianCS.C

\*---------------------------------------------------------------------------*/

#ifndef cartesianCS_H
#define cartesianCS_H

#include "coordinateSystem.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class cartesianCS Declaration
\*---------------------------------------------------------------------------*/

class cartesianCS
:
    public coordinateSystem
{
protected:

    // Protected Member Functions


        //- Convert from local coordinate system to the global Cartesian system
        //  with optional translation for the origin
        virtual vector localToGlobal(const vector&, bool translate) const;

        //- Convert from local coordinate system to the global Cartesian system
        //  with optional translation for the origin
        virtual tmp<vectorField> localToGlobal
        (
            const vectorField&,
            bool translate
        ) const;

        //- Convert from global Cartesian system to the local coordinate system
        //  with optional translation for the origin
        virtual vector globalToLocal(const vector&, bool translate) const;

        //- Convert from global Cartesian system to the local coordinate system
        //  with optional translation for the origin
        virtual tmp<vectorField> globalToLocal
        (
            const vectorField&,
            bool translate
        ) const;


public:

    //- Runtime type information
    TypeName("cartesian");


    // Constructors

        //- Construct null
        cartesianCS();

        //- Construct copy
        cartesianCS
        (
            const coordinateSystem&
        );

        //- Construct copy with a different name
        cartesianCS
        (
            const word& name,
            const coordinateSystem&
        );

        //- Construct from origin and rotation
        cartesianCS
        (
            const word& name,
            const point& origin,
            const coordinateRotation&
        );

        //- Construct from origin and 2 axes
        cartesianCS
        (
            const word& name,
            const point& origin,
            const vector& axis,
            const vector& dirn
        );

        //- Construct from dictionary
        cartesianCS(const word&, const dictionary&);


        //- Construct from dictionary and objectRegistry
        cartesianCS(const objectRegistry&, const dictionary&);


    //- Destructor
    virtual ~cartesianCS();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
