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
    Foam::FunctionalDimensionedField

Description
    DimensionedField with a corresponding run-time selected function
    to evaluate and update the field.

See also
    Foam::DimensionedField
    Foam::DimensionedFieldFunction
    Foam::UFunctionalDimensionedField

SourceFiles
    FunctionalDimensionedField.C

\*---------------------------------------------------------------------------*/

#ifndef FunctionalDimensionedField_H
#define FunctionalDimensionedField_H

#include "FunctionalDimensionedFieldFwd.H"
#include "DimensionedField.H"
#include "DimensionedFieldFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions
template<class Type, class GeoMesh>
void writeEntry
(
    Ostream& os,
    const FunctionalDimensionedField<Type, GeoMesh>& udff
);

/*---------------------------------------------------------------------------*\
                    Class FunctionalDimensionedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
class FunctionalDimensionedField
:
    public DimensionedField<Type, GeoMesh>
{
    // Private Data

        //- Name of the function entry/sub dictionary
        const word funcName_;

        //- Optional function pointer
        autoPtr<DimensionedFieldFunction<DimensionedField<Type, GeoMesh>>>
            funcPtr_;


public:

    // Constructors

        //- Construct from name, mesh, dimensions, field and dictionary
        FunctionalDimensionedField
        (
            const word& name,
            const word& funcName,
            const GeoMesh&,
            const dimensionSet& dimensions,
            const dictionary&
        );

        //- Construct from name, mesh, dimensions, field, dictionary
        //  and default value
        FunctionalDimensionedField
        (
            const word& name,
            const word& funcName,
            const GeoMesh&,
            const dimensionSet& dimensions,
            const dictionary&,
            const Type& defaultValue
        );

        //- Construct from given FunctionalDimensionedField with new mesh
        FunctionalDimensionedField
        (
            const FunctionalDimensionedField<Type, GeoMesh>&,
            const GeoMesh&
        );

        //- Construct by mapping given FunctionalDimensionedField
        FunctionalDimensionedField
        (
            const FunctionalDimensionedField<Type, GeoMesh>&,
            const GeoMesh&,
            const fieldMapper&
        );

        //- Copy constructor
        FunctionalDimensionedField
        (
            const FunctionalDimensionedField<Type, GeoMesh>&
        );


    // Member Functions

        //- Map if mapping direct or evaluate if function defined
        void map
        (
            const DimensionedField<Type, GeoMesh>& df,
            const fieldMapper& mapper
        );

        //- Reset
        void reset();

        //- Update
        bool update();

        using regIOobject::write;

        //- Write
        void write(Ostream&) const;


    // Global functions

        friend void writeEntry <Type, GeoMesh>
        (
            Ostream& os,
            const FunctionalDimensionedField<Type, GeoMesh>& udff
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FunctionalDimensionedField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
