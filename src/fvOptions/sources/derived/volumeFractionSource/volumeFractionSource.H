/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::fv::volumeFractionSource

Description
    This option adds transport terms into the equations to account for the
    presence of a constant volume fraction. The volume fraction is read from
    constant/alpha.<phase>, where <phase> is given as a parameter to the
    option. Both advective and diffusive terms are added, and the resulting
    solution is time-accurate. The flux and velocity are treated as
    superficial.

    This can be used to represent the effect of porous media that are caused
    purely by the reduction in volume of the fluid phase; i.e., additional
    blockage, and changes to transport and diffusion rates. It does not
    represent losses or transfers with the porous media. That requires separate
    sub-modelling.

Usage
    \table
        Property    | Description                            | Req'd? | Default
        phase       | Name of the phase associated with the volume fraction \\
                                                             | yes    |
        phi         | Name of the flux field                 | no     | phi
        rho         | Name of the density field              | no     | rho
        U           | Name of the velocity field             | no     | U
        fields      | Names of the fields to apply the option to \\
                                                             | yes    |
    \endtable

    Example specification:
    \verbatim
    <fvOptionName>
    {
        type        volumeFractionSource;
        phase       solid;
        phi         phi;
        rho         rho;
        U           U;
        fields      (rho U e);
    }
    \endverbatim

SourceFiles
    volumeFractionSource.C

\*---------------------------------------------------------------------------*/

#ifndef volumeFractionSource_H
#define volumeFractionSource_H

#include "fvOption.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class volumeFractionSource Declaration
\*---------------------------------------------------------------------------*/

class volumeFractionSource
:
    public option
{
    // Private Member Data

        //- The name of the phase
        const word phaseName_;

        //- The name of the flux field
        word phiName_;

        //- The name of the density field
        word rhoName_;

        //- The name of the velocity field
        word UName_;


    // Private Member Functions

        //- Get the volume fraction field
        const volScalarField& alpha() const;

        //- Get the diffusivity for a given field
        tmp<volScalarField> D(const label) const;

        //- Add time-derivative terms to an equation
        template <class Type>
        void addDdtSup(fvMatrix<Type>&, const label);

        //- Add time-derivative terms to a compressible equation
        template <class Type>
        void addDdtSup(const volScalarField&, fvMatrix<Type>&, const label);

        //- Add time-derivative terms to a phase-compressible equation
        template <class Type>
        void addDdtSup
        (
            const volScalarField&,
            const volScalarField&,
            fvMatrix<Type>&,
            const label
        );

        //- Add divergence terms to an equation
        template <class Type>
        void addDivSup(fvMatrix<Type>&, const label);

        //- Add divergence terms to the momentum equation
        void addUDivSup(fvMatrix<vector>&, const label);

        //- Add divergence terms to the continuity equation
        void addRhoDivSup(fvMatrix<scalar>&, const label);

        //- Add Laplacian terms to an equation
        template <class Type, class AlphaFieldType>
        void addLaplacianSup
        (
            const AlphaFieldType& alpha,
            fvMatrix<Type>& eqn,
            const label fieldi
        );


public:

    //- Runtime type information
    TypeName("volumeFractionSource");


    // Constructors

        //- Construct from components
        volumeFractionSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        volumeFractionSource(const volumeFractionSource&) = delete;


    //- Destructor
    virtual ~volumeFractionSource();


    // Member Functions

        // Evaluation

            // Explicit and implicit sources

                virtual void addSup(fvMatrix<scalar>&, const label);

                virtual void addSup(fvMatrix<vector>&, const label);

                virtual void addSup(fvMatrix<symmTensor>&, const label);

                virtual void addSup(fvMatrix<sphericalTensor>&, const label);

                virtual void addSup(fvMatrix<tensor>&, const label);


            // Explicit and implicit sources for compressible equations

                virtual void addSup
                (
                    const volScalarField&,
                    fvMatrix<scalar>&,
                    const label
                );

                virtual void addSup
                (
                    const volScalarField&,
                    fvMatrix<vector>&,
                    const label
                );

                virtual void addSup
                (
                    const volScalarField&,
                    fvMatrix<symmTensor>&,
                    const label
                );

                virtual void addSup
                (
                    const volScalarField&,
                    fvMatrix<sphericalTensor>&,
                    const label
                );

                virtual void addSup
                (
                    const volScalarField&,
                    fvMatrix<tensor>&,
                    const label
                );


            // Explicit and implicit sources for phase equations

                virtual void addSup
                (
                    const volScalarField&,
                    const volScalarField&,
                    fvMatrix<scalar>&,
                    const label
                );

                virtual void addSup
                (
                    const volScalarField&,
                    const volScalarField&,
                    fvMatrix<vector>&,
                    const label
                );

                virtual void addSup
                (
                    const volScalarField&,
                    const volScalarField&,
                    fvMatrix<symmTensor>&,
                    const label
                );

                virtual void addSup
                (
                    const volScalarField&,
                    const volScalarField&,
                    fvMatrix<sphericalTensor>&,
                    const label
                );

                virtual void addSup
                (
                    const volScalarField&,
                    const volScalarField&,
                    fvMatrix<tensor>&,
                    const label
                );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const volumeFractionSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
