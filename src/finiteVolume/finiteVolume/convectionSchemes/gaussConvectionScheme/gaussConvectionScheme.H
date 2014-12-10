/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::fv::gaussConvectionScheme

Description
    Basic second-order convection using face-gradients and Gauss' theorem.

SourceFiles
    gaussConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef gaussConvectionScheme_H
#define gaussConvectionScheme_H

#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

//- Temporary debug switch to provide warning about backward-compatibility
//  issue with setting div schemes for steady-state
extern int warnUnboundedGauss;

/*---------------------------------------------------------------------------*\
                       Class gaussConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gaussConvectionScheme
:
    public fv::convectionScheme<Type>
{
    // Private data

        tmp<surfaceInterpolationScheme<Type> > tinterpScheme_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        gaussConvectionScheme(const gaussConvectionScheme&);

        //- Disallow default bitwise assignment
        void operator=(const gaussConvectionScheme&);


public:

    //- Runtime type information
    TypeName("Gauss");


    // Constructors

        //- Construct from flux and interpolation scheme
        gaussConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const tmp<surfaceInterpolationScheme<Type> >& scheme
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            tinterpScheme_(scheme)
        {}

        //- Construct from flux and Istream
        gaussConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            tinterpScheme_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {
            is.rewind();
            word bounded(is);

            if
            (
                warnUnboundedGauss
             && word(mesh.ddtScheme("default")) == "steadyState"
             && bounded != "bounded"
            )
            {
                fileNameList controlDictFiles(findEtcFiles("controlDict"));

                IOWarningIn("gaussConvectionScheme", is)
                    << "Unbounded 'Gauss' div scheme used in "
                       "steady-state solver, use 'bounded Gauss' "
                       "to ensure boundedness.\n"
                    << "    To remove this warning switch off "
                    << "'warnUnboundedGauss' in "
                    << controlDictFiles[controlDictFiles.size()-1]
                    << endl;
            }
        }


    // Member Functions

        const surfaceInterpolationScheme<Type>& interpScheme() const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > flux
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<fvMatrix<Type> > fvmDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "gaussConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
