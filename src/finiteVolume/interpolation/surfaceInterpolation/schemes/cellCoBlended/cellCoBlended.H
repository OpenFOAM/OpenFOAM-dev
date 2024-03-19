/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2024 OpenFOAM Foundation
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
    Foam::cellCoBlended

Description
    Two-scheme cell-based Courant number based blending interpolation scheme.

    This scheme is equivalent to the CoBlended scheme except that the Courant
    number is evaluated for cells using the same approach as use in the
    finite-volume solvers and then interpolated to the faces rather than being
    estimated directly at the faces based on the flux.  This is a more
    consistent method for evaluating the Courant number but suffers from the
    need to interpolate which introduces a degree of freedom.  However, the
    interpolation scheme for "Co" is run-time selected and may be specified in
    "interpolationSchemes" and "localMax" might be most appropriate.

    Example of the cellCoBlended scheme specification using LUST for Courant
    numbers less than 1 and linearUpwind for Courant numbers greater than 10:
    \verbatim
    divSchemes
    {
        .
        .
        div(phi,U)  Gauss cellCoBlended 1 LUST grad(U) 10 linearUpwind grad(U);
        .
        .
    }

    interpolationSchemes
    {
        .
        .
        interpolate(Co) localMax;
        .
        .
    }
    \endverbatim

See also
    Foam::CoBlended
    Foam::localBlended

SourceFiles
    cellCoBlended.C

\*---------------------------------------------------------------------------*/

#ifndef cellCoBlended_H
#define cellCoBlended_H

#include "surfaceInterpolationScheme.H"
#include "blendedSchemeBase.H"
#include "surfaceInterpolate.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cellCoBlended Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cellCoBlended
:
    public surfaceInterpolationScheme<Type>,
    public blendedSchemeBase<Type>
{
    // Private Data

        //- Courant number below which scheme1 is used
        const scalar Co1_;

        //- Scheme 1
        tmp<surfaceInterpolationScheme<Type>> tScheme1_;

        //- Courant number above which scheme2 is used
        const scalar Co2_;

        //- Scheme 2
        tmp<surfaceInterpolationScheme<Type>> tScheme2_;

        //- The face-flux used to compute the face Courant number
        const surfaceScalarField& faceFlux_;


public:

    //- Runtime type information
    TypeName("cellCoBlended");


    // Constructors

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        cellCoBlended
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            Co1_(readScalar(is)),
            tScheme1_
            (
                surfaceInterpolationScheme<Type>::New(mesh, is)
            ),
            Co2_(readScalar(is)),
            tScheme2_
            (
                surfaceInterpolationScheme<Type>::New(mesh, is)
            ),
            faceFlux_
            (
                mesh.lookupObject<surfaceScalarField>(word(is))
            )
        {
            if (Co1_ < 0 || Co2_ < 0 || Co1_ >= Co2_)
            {
                FatalIOErrorInFunction(is)
                    << "coefficients = " << Co1_ << " and " << Co2_
                    << " should be > 0 and Co2 > Co1"
                    << exit(FatalIOError);
            }
        }


        //- Construct from mesh, faceFlux and Istream
        cellCoBlended
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            Co1_(readScalar(is)),
            tScheme1_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            ),
            Co2_(readScalar(is)),
            tScheme2_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            ),
            faceFlux_(faceFlux)
        {
            if (Co1_ < 0 || Co2_ < 0 || Co1_ >= Co2_)
            {
                FatalIOErrorInFunction(is)
                    << "coefficients = " << Co1_ << " and " << Co2_
                    << " should be > 0 and Co2 > Co1"
                    << exit(FatalIOError);
            }
        }

        //- Disallow default bitwise copy construction
        cellCoBlended(const cellCoBlended&) = delete;


    // Member Functions

        //- Return the face-based blending factor
        virtual tmp<surfaceScalarField> blendingFactor
        (
             const VolField<Type>& vf
        ) const
        {
            const fvMesh& mesh = this->mesh();
            tmp<surfaceScalarField> tUflux = faceFlux_;

            if (faceFlux_.dimensions() == dimMassFlux)
            {
                // Currently assume that the density field
                // corresponding to the mass-flux is named "rho"
                const volScalarField& rho =
                    mesh.objectRegistry::template lookupObject<volScalarField>
                    ("rho");

                tUflux = faceFlux_/fvc::interpolate(rho);
            }
            else if (faceFlux_.dimensions() != dimFlux)
            {
                FatalErrorInFunction
                    << "dimensions of faceFlux are not correct"
                    << exit(FatalError);
            }

            volScalarField Co
            (
                IOobject
                (
                    "Co",
                    mesh.time().name(),
                    mesh
                ),
                mesh,
                dimensionedScalar(dimless, 0),
                extrapolatedCalculatedFvPatchScalarField::typeName
            );

            scalarField sumPhi
            (
                fvc::surfaceSum(mag(tUflux))().primitiveField()
            );

            Co.primitiveFieldRef() =
                (sumPhi/mesh.V().primitiveField())
               *(0.5*mesh.time().deltaTValue());

            Co.correctBoundaryConditions();

            return surfaceScalarField::New
            (
                vf.name() + "BlendingFactor",
                scalar(1)
              - max
                (
                    min
                    (
                        (fvc::interpolate(Co) - Co1_)/(Co2_ - Co1_),
                        scalar(1)
                    ),
                    scalar(0)
                )
            );
        }


        //- Return the interpolation weighting factors
        tmp<surfaceScalarField>
        weights
        (
            const VolField<Type>& vf
        ) const
        {
            surfaceScalarField bf(blendingFactor(vf));

            return
                bf*tScheme1_().weights(vf)
              + (scalar(1) - bf)*tScheme2_().weights(vf);
        }


        //- Return the face-interpolate of the given cell field
        //  with explicit correction
        tmp<SurfaceField<Type>>
        interpolate
        (
            const VolField<Type>& vf
        ) const
        {
            surfaceScalarField bf(blendingFactor(vf));

            return
                bf*tScheme1_().interpolate(vf)
              + (scalar(1) - bf)*tScheme2_().interpolate(vf);
        }


        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return tScheme1_().corrected() || tScheme2_().corrected();
        }


        //- Return the explicit correction to the face-interpolate
        //  for the given field
        virtual tmp<SurfaceField<Type>>
        correction
        (
            const VolField<Type>& vf
        ) const
        {
            surfaceScalarField bf(blendingFactor(vf));

            if (tScheme1_().corrected())
            {
                if (tScheme2_().corrected())
                {
                    return
                    (
                        bf
                      * tScheme1_().correction(vf)
                      + (scalar(1) - bf)
                      * tScheme2_().correction(vf)
                    );
                }
                else
                {
                    return
                    (
                        bf
                      * tScheme1_().correction(vf)
                    );
                }
            }
            else if (tScheme2_().corrected())
            {
                return
                (
                    (scalar(1) - bf)
                  * tScheme2_().correction(vf)
                );
            }
            else
            {
                return tmp<SurfaceField<Type>>
                (
                    nullptr
                );
            }
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cellCoBlended&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
