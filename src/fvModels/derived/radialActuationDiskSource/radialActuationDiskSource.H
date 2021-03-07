/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::fv::radialActuationDiskSource

Description
    Actuation disk source including radial thrust

    Constant values for momentum source for actuation disk

        \f[
            T = 2 \rho A U_{o}^2 a (1-a)
        \f]

    and

        \f[
            U_1 = (1 - a)U_{o}
        \f]

    where:
    \vartable
        A   | disk area
        U_o | upstream velocity
        a   | 1 - Cp/Ct
        U_1 | velocity at the disk
    \endvartable

    The thrust is distributed by a radial function:

        \f[
            thrust(r) = T (C_0 + C_1 r^2 + C_2 r^4)
        \f]

Usage
    Example usage:
    \verbatim
    radialActuationDiskSource1
    {
        type            radialActuationDiskSource;

        selectionMode   cellSet;
        cellSet         radialActuationDisk1;

        U               U;              // Name of the velocity field

        diskDir         (-1 0 0);       // Disk direction
        Cp              0.1;            // Power coefficient
        Ct              0.5;            // Thrust coefficient
        diskArea        5.0;            // Disk area
        coeffs          (0.1 0.5 0.01); // Radial distribution coefficients
        upstreamPoint   (0 0 0);        // Upstream point
    }
    \endverbatim


SourceFiles
    radialActuationDiskSource.C
    radialActuationDiskSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef radialActuationDiskSource_H
#define radialActuationDiskSource_H

#include "actuationDiskSource.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class radialActuationDiskSource Declaration
\*---------------------------------------------------------------------------*/

class radialActuationDiskSource
:
    public actuationDiskSource
{
    // Private Data

        //- Coefficients for the radial distribution
        FixedList<scalar, 3> radialCoeffs_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- Add resistance to the UEqn
        template<class RhoFieldType>
        void addRadialActuationDiskAxialInertialResistance
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;


public:

    //- Runtime type information
    TypeName("radialActuationDiskSource");


    // Constructors

        //- Construct from components
        radialActuationDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        radialActuationDiskSource(const radialActuationDiskSource&) = delete;


    //- Destructor
    virtual ~radialActuationDiskSource()
    {}


    // Member Functions

        //- Source term to momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const word& fieldName
        ) const;

        //- Source term to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const word& fieldName
        ) const;


    // IO

        //- Read dictionary
        virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const radialActuationDiskSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "radialActuationDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
