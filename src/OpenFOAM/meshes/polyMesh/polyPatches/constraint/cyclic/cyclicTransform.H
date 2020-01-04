/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::cyclicTransform

Description
    Cyclic plane tranformation.

SourceFiles
    cyclicTransform.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicTransform_H
#define cyclicTransform_H

#include "transformer.H"
#include "dictionary.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cyclicTransform Declaration
\*---------------------------------------------------------------------------*/

class cyclicTransform
{
public:

    enum transformTypes
    {
        UNSPECIFIED,        // Unspecified -> automatic transformation
        NONE,               // No tranformation
        ROTATIONAL,         // Rotation around a coordinate axis
        TRANSLATIONAL       // Translation
    };

    static const NamedEnum<transformTypes, 4> transformTypeNames;


private:

        //- Type of transformation
        transformTypes transformType_;


protected:

        //- Transformation between the coupled patches
        mutable transformer transform_;


public:

    // Constructors

        // Unspecified tranformation
        cyclicTransform();

        // Transformation read from dictionary
        cyclicTransform(const dictionary& dict);


    // Member Functions

        // Access

            //- Type of transform
            transformTypes transformType() const
            {
                return transformType_;
            }

            //- Return transformation between the coupled patches
            const transformer& transform() const
            {
                return transform_;
            }


        //- Write the polyPatch data as a dictionary
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
