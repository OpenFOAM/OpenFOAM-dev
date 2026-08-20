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

Class
    Foam::DimensionedFieldFunctions::Surfaces

Description
    Dimensioned field function which sets constant values inside a set of
    surfaces, and a default value elsewhere. Cutting is used to set fractional
    values in elements that intersect the surfaces. Foam::searchableSurface is
    used to define the surface geometry, so triangulated surface files can be
    used, as well as a list of simple predefined shapes (box, cylinder, ...).

    Example specification:
    internalField
    {
        type        surfaces;

        // Set a value of one everywhere, except ...
        defaultValue 1;

        surfaces
        {
            // Set a value of zero in a sphere
            bubble
            {
                type        sphere;

                centre      (0.5 0.5 0.5);
                radius      0.1;

                value       0;
            }

            // Set a value of zero where the y coordinate is greater than 1
            air
            {
                type        plane;

                point       (0 1 0);
                normal      (0 -1 0);

                value       0;
            }
        }
    }

SourceFiles
    Surfaces_DimensionedFieldFunction.C

\*---------------------------------------------------------------------------*/

#ifndef Surfaces_DimensionedFieldFunction_H
#define Surfaces_DimensionedFieldFunction_H

#include "DimensionedFieldFunction.H"
#include "searchableSurfaceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class polyPatch;
class fvMesh;
class fvPatch;

namespace DimensionedFieldFunctions
{

/*---------------------------------------------------------------------------*\
                            Class Surfaces Declaration
\*---------------------------------------------------------------------------*/

template<class DimensionedFieldType>
class Surfaces
:
    public DimensionedFieldFunction<DimensionedFieldType>
{
private:

    // Private Type Definitions

        //- The primitive data type
        typedef typename DimensionedFieldType::Type_ Type;


    // Private Data

        //- Default (outside) value
        dimensioned<Type> defaultValue_;

        //- Surfaces dictionary
        dictionary surfacesDict_;

        //- Surfaces
        autoPtr<searchableSurfaceList> surfacesPtr_;

        //- Inside values
        PtrList<dimensioned<Type>> values_;


    // Private Member Functions

        //- Return polyMesh given an fvMesh
        static const polyMesh& poly(const fvMesh&);

        //- Return polyPatch given an fvPatch
        static const polyPatch& poly(const fvPatch&);

        //- Read the surfaces and the values
        void readSurfacesAndValues();


public:

    //- Runtime type information
    TypeName("surfaces");


    // Constructors

        //- Construct with dictionary to initialise given field
        Surfaces
        (
            const dictionary& dict,
            DimensionedFieldType& field
        );

        //- Construct a copy for the given field
        Surfaces
        (
            const Surfaces& dff,
            DimensionedFieldType& field
        );

        //- Construct and return a clone for the specified field
        virtual autoPtr<DimensionedFieldFunction<DimensionedFieldType>> clone
        (
            DimensionedFieldType& field
        ) const;


    //- Destructor
    virtual ~Surfaces()
    {}


    // Member Functions

        //- Evaluate the function and set the field
        virtual void evaluate();

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace DimensionedFieldFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Surfaces_DimensionedFieldFunction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
