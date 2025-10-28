/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::fv::bXiKernelCorr

Description
    Early kernel propagation rate correction based on specified shape

Usage
    Example usage:
    \verbatim
    bXiKernelCorr1
    {
        type            bXiKernelCorr;

        cellZone        kernelCorrection;

        ignition        ignition;

        duration        0.003;

        kernelShape
        {
            type            spherical;
        }
    }
    \endverbatim

    Where:
    \table
        Property    | Description                 | Required   | Default value
        cellZone    | Correction cellZone         | yes        |
        ignition    | Name of the ignition fvModel to correct | yes |
        duration    | Maximum correction duration | yes        |
        kernelShape | Kernel shape                | yes        |
        bMin        | Min b above which the correction is applied | no | \\
                      XiFluid::bMin/2
    \endtable

    The correction is applied if \f$ min(b) > bMin \f$ and time of application
    is less than the specified \c duration.

See also
    Foam::kernelShape
    Foam::kernelShapes::planar
    Foam::kernelShapes::cylindrical
    Foam::kernelShapes::spherical

SourceFiles
    bXiKernelCorr.C

\*---------------------------------------------------------------------------*/

#ifndef bXiKernelCorr_H
#define bXiKernelCorr_H

#include "bXiIgnition.H"
#include "kernelShape.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class bXiKernelCorr Declaration
\*---------------------------------------------------------------------------*/

class bXiKernelCorr
:
    public fvModel
{

protected:

    // Protected Data

        //- The correction cellZone
        fvCellZone zone_;

        autoPtr<kernelShape> kernelShape_;

        const fv::bXiIgnition& ignition_;

        //- Maximum correction duration
        dimensionedScalar duration_;

        //- Minimum b below which the Xi correction is not applied
        //  Must be less that XiFluid::bMin and defaults to XiFluid::bMin/2
        dimensionedScalar bMin_;


private:

    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("bXiKernelCorr");


    // Constructors

        //- Construct from explicit source name and mesh
        bXiKernelCorr
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        bXiKernelCorr
        (
            const bXiKernelCorr&
        ) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Add explicit and implicit contributions to compressible equation

            //- Add ignition contribution to b equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& b,
                fvMatrix<scalar>& eqn
            ) const;


        // Mesh motion

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const bXiKernelCorr&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
