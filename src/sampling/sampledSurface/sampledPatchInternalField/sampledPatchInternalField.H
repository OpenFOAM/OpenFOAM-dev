/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::sampledSurfaces::patchInternalField

Description
    Variation of sampledPatch that samples the internalField (at a given
    normal distance from the patch) instead of the patchField.
    Note:
    - interpolate=false : get cell value on faces
    - interpolate=true  : interpolate inside cell and interpolate to points
    There is no option to get interpolated value inside the cell on the faces.

SourceFiles
    sampledPatchInternalField.C

\*---------------------------------------------------------------------------*/

#ifndef sampledPatchInternalField_H
#define sampledPatchInternalField_H

#include "sampledPatch.H"
#include "mappedInternalPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSurfaces
{

/*---------------------------------------------------------------------------*\
                     Class patchInternalField Declaration
\*---------------------------------------------------------------------------*/

class patchInternalField
:
    public patch
{
    // Private Data

        //- Mapping engines
        PtrList<mappedInternalPatchBase> mappers_;


    // Private Member Functions

        //- Sample field on faces
        template<class Type>
        tmp<Field<Type>> sampleField
        (
            const VolField<Type>& vField
        ) const;

        template<class Type>
        tmp<Field<Type>> interpolateField(const interpolation<Type>&) const;


public:

    //- Runtime type information
    TypeName("patchInternalField");


    // Constructors

        //- Construct from dictionary
        patchInternalField
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~patchInternalField();


    // Member Functions

        // Sample

            //- Sample field on surface
            virtual tmp<scalarField> sample
            (
                const volScalarField&
            ) const;

            //- Sample field on surface
            virtual tmp<vectorField> sample
            (
                const volVectorField&
            ) const;

            //- Sample field on surface
            virtual tmp<sphericalTensorField> sample
            (
                const volSphericalTensorField&
            ) const;

            //- Sample field on surface
            virtual tmp<symmTensorField> sample
            (
                const volSymmTensorField&
            ) const;

            //- Sample field on surface
            virtual tmp<tensorField> sample
            (
                const volTensorField&
            ) const;


        // Interpolate

            //- Interpolate field on surface
            virtual tmp<scalarField> interpolate
            (
                const interpolation<scalar>&
            ) const;

            //- Interpolate field on surface
            virtual tmp<vectorField> interpolate
            (
                const interpolation<vector>&
            ) const;

            //- Interpolate field on surface
            virtual tmp<sphericalTensorField> interpolate
            (
                const interpolation<sphericalTensor>&
            ) const;

            //- Interpolate field on surface
            virtual tmp<symmTensorField> interpolate
            (
                const interpolation<symmTensor>&
            ) const;

            //- Interpolate field on surface
            virtual tmp<tensorField> interpolate
            (
                const interpolation<tensor>&
            ) const;


        //- Write
        virtual void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sampledSurfaces
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledPatchInternalFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
