/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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
    Foam::nonConformalMappedPolyFacesFvsPatchLabelField

SourceFiles
    nonConformalMappedPolyFacesFvsPatchLabelField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalMappedPolyFacesFvsPatchLabelField_H
#define nonConformalMappedPolyFacesFvsPatchLabelField_H

#include "fvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class nonConformalMappedPolyFacesFvsPatchLabel Declaration
\*---------------------------------------------------------------------------*/

class nonConformalMappedPolyFacesFvsPatchLabelField
:
    public fvsPatchLabelField
{
private:

    // Private Member Data

        //- The offsets of the processor sub-ranges
        labelList procOffsets_;


public:

    //- Runtime type information
    TypeName("nonConformalMappedPolyFaces");


    // Constructors

        //- Construct from patch and internal field
        nonConformalMappedPolyFacesFvsPatchLabelField
        (
            const fvPatch&,
            const DimensionedField<label, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nonConformalMappedPolyFacesFvsPatchLabelField
        (
            const fvPatch&,
            const DimensionedField<label, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        nonConformalMappedPolyFacesFvsPatchLabelField
        (
            const nonConformalMappedPolyFacesFvsPatchLabelField&,
            const fvPatch&,
            const DimensionedField<label, surfaceMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        nonConformalMappedPolyFacesFvsPatchLabelField
        (
            const nonConformalMappedPolyFacesFvsPatchLabelField&
        ) = delete;

        //- Copy constructor setting internal field reference
        nonConformalMappedPolyFacesFvsPatchLabelField
        (
            const nonConformalMappedPolyFacesFvsPatchLabelField&,
            const DimensionedField<label, surfaceMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvsPatchLabelField> clone
        (
            const DimensionedField<label, surfaceMesh>& iF
        ) const
        {
            return tmp<fvsPatchLabelField>
            (
                new nonConformalMappedPolyFacesFvsPatchLabelField(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvsPatchLabelField&);


        //- Return the offsets of the processor sub-ranges
        inline const labelList& procOffsets() const
        {
            return procOffsets_;
        }

        //- Access the offsets of the processor sub-ranges
        inline labelList& procOffsets()
        {
            return procOffsets_;
        }

        //- Return the sizes of the processor sub-ranges
        labelList procSizes() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
