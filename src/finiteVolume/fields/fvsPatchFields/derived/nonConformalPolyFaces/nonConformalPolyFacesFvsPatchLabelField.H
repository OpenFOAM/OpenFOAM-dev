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
    Foam::nonConformalPolyFacesFvsPatchLabelField

SourceFiles
    nonConformalPolyFacesFvsPatchLabelField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalPolyFacesFvsPatchLabelField_H
#define nonConformalPolyFacesFvsPatchLabelField_H

#include "fvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class nonConformalPolyFacesFvsPatchLabel Declaration
\*---------------------------------------------------------------------------*/

class nonConformalPolyFacesFvsPatchLabelField
:
    public fvsPatchLabelField
{

public:

    //- Runtime type information
    TypeName("nonConformalPolyFaces");


    // Constructors

        //- Construct from patch and internal field
        nonConformalPolyFacesFvsPatchLabelField
        (
            const fvPatch&,
            const DimensionedField<label, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nonConformalPolyFacesFvsPatchLabelField
        (
            const fvPatch&,
            const DimensionedField<label, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        nonConformalPolyFacesFvsPatchLabelField
        (
            const nonConformalPolyFacesFvsPatchLabelField&,
            const fvPatch&,
            const DimensionedField<label, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        nonConformalPolyFacesFvsPatchLabelField
        (
            const nonConformalPolyFacesFvsPatchLabelField&
        ) = delete;

        //- Copy constructor setting internal field reference
        nonConformalPolyFacesFvsPatchLabelField
        (
            const nonConformalPolyFacesFvsPatchLabelField&,
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
                new nonConformalPolyFacesFvsPatchLabelField(*this, iF)
            );
        }


    // Member Functions

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
