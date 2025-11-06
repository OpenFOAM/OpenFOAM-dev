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
    Foam::fv::GaussianbXiIgnition

Description
    Gaussian profile time dependent rate ignition model for the Weller b-Xi
    combustion model

Usage
    Example usage:
    \verbatim
    GaussianbXiIgnition1
    {
        type            GaussianbXiIgnition;

        cellZone        ignition;

        position        (0 0 0);
        diameter        1e-3;

        kernelShape
        {
            type            spherical;
        }

        start           0;
        duration        0.003;
        rate            1e3;
    }
    \endverbatim

    Where:
    \table
        Property    | Description               | Required   | Default value
        cellZone    | Correction cellZone       | yes        |
        position    | Ignition centre           | yes        |
        diameter    | Diameter of the ignition kernel | yes        |
        kernelShape | Kernel shape              | yes        |
        start       | Ignition start time       | yes        |
        duration    | Ignition duration         | yes        |
        rate        | Burning rate coefficient Function1 of time [1/s] | yes |
    \endtable

See also
    Foam::fv::constantbXiIgnition
    Foam::fv::bXiKernelCorr
    Foam::kernelShape
    Foam::kernelShapes::planar
    Foam::kernelShapes::cylindrical
    Foam::kernelShapes::spherical

SourceFiles
    GaussianbXiIgnition.C

\*---------------------------------------------------------------------------*/

#ifndef GaussianbXiIgnition_H
#define GaussianbXiIgnition_H

#include "bXiTimedIgnition.H"
#include "kernelShape.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class GaussianbXiIgnition Declaration
\*---------------------------------------------------------------------------*/

class GaussianbXiIgnition
:
    public bXiTimedIgnition
{

protected:

    // Protected Data

        //- The ignition cellZone
        fvCellZone zone_;

        autoPtr<kernelShape> kernelShape_;

        dimensionedVector position_;

        dimensionedScalar diameter_;

        //- Ignition burning rate [1/t]
        autoPtr<Function1<scalar>> rate_;

        //- Kernel volume correction factor
        scalarList Vk_;


private:

    // Private Member Functions

        //- Calculate the kernel volume correction factor
        void calcVk();

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("GaussianbXiIgnition");


    // Constructors

        //- Construct from explicit source name and mesh
        GaussianbXiIgnition
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        GaussianbXiIgnition
        (
            const GaussianbXiIgnition&
        ) = delete;


    // Member Functions

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
        void operator=(const GaussianbXiIgnition&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
