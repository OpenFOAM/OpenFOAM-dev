/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::polyFacesFvsPatchLabelField

SourceFiles
    polyFacesFvsPatchLabelField.C

\*---------------------------------------------------------------------------*/

#ifndef polyFacesFvsPatchLabelField_H
#define polyFacesFvsPatchLabelField_H

#include "fvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class polyFacesFvsPatchLabelField Declaration
\*---------------------------------------------------------------------------*/

class polyFacesFvsPatchLabelField
:
    public fvsPatchLabelField
{
private:

    // Private Member Functions

        //- Initialise
        void init();


public:

    //- Runtime type information
    TypeName("polyFaces");


    // Constructors

        //- Construct from patch and internal field
        polyFacesFvsPatchLabelField
        (
            const fvPatch&,
            const DimensionedField<label, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        polyFacesFvsPatchLabelField
        (
            const fvPatch&,
            const DimensionedField<label, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        polyFacesFvsPatchLabelField
        (
            const polyFacesFvsPatchLabelField&,
            const fvPatch&,
            const DimensionedField<label, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        polyFacesFvsPatchLabelField
        (
            const polyFacesFvsPatchLabelField&
        ) = delete;

        //- Copy constructor setting internal field reference
        polyFacesFvsPatchLabelField
        (
            const polyFacesFvsPatchLabelField&,
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
                new polyFacesFvsPatchLabelField(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
