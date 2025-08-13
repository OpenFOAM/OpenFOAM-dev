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
    Foam::totalPressureVelocityMagnitudeLagrangianScalarFieldSource

Description
    Mix-in for source conditions that compute the velocity magnitude from the
    difference between a specified upstream pressure and the local interpolated
    pressure

SourceFiles
    totalPressureVelocityMagnitudeLagrangianScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef totalPressureVelocityMagnitudeLagrangianScalarFieldSource_H
#define totalPressureVelocityMagnitudeLagrangianScalarFieldSource_H

#include "LagrangianInjection.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class LagrangianFieldSourceBase;
class cloudLagrangianFieldSource;

/*---------------------------------------------------------------------------*\
  Class totalPressureVelocityMagnitudeLagrangianScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class totalPressureVelocityMagnitudeLagrangianScalarFieldSource
{
private:

    // Private Data

        //- Reference to the field source
        const LagrangianFieldSourceBase& field_;

        //- Reference to the cloud field source
        const cloudLagrangianFieldSource& cloudField_;

        //- Total pressure entry. Cached until evaluation as the function
        //  cannot be constructed until the carrier pressure dimensions are
        //  available. After this it is deleted.
        mutable autoPtr<entry> p0Entry_;

        //- Total pressure
        mutable autoPtr<Function1<scalar>> p0_;

        //- Name of the carrier density field
        const word rhocName_;

        //- Name of the carrier pressure field
        const word pcName_;


public:

    // Constructors

        //- Construct from a dictionary
        totalPressureVelocityMagnitudeLagrangianScalarFieldSource
        (
            const LagrangianFieldSourceBase&,
            const cloudLagrangianFieldSource&,
            const dictionary& dict
        );

        //- Copy constructor
        totalPressureVelocityMagnitudeLagrangianScalarFieldSource
        (
            const totalPressureVelocityMagnitudeLagrangianScalarFieldSource&,
            const LagrangianFieldSourceBase&,
            const cloudLagrangianFieldSource&
        );


    //- Destructor
    ~totalPressureVelocityMagnitudeLagrangianScalarFieldSource();


    // Member Functions

        //- Return the value for an instantaneous injection
        tmp<LagrangianSubScalarField> Umag
        (
            const LagrangianInjection&,
            const LagrangianSubMesh&
        ) const;

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
