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
    Foam::uniformSizeNumberLagrangianScalarFieldSource

Description
    Base class for Lagrangian source conditions that calculate the number field
    from a total (e.g., mass, volume or flowRate) and a uniform size of the
    created parcels (e.g., mass, volume, surface area or number).

SourceFiles
    uniformSizeNumberLagrangianScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef uniformSizeNumberLagrangianScalarFieldSource_H
#define uniformSizeNumberLagrangianScalarFieldSource_H

#include "LagrangianFieldSources.H"
#include "CloudLagrangianFieldSource.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class uniformSizeNumberLagrangianScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class uniformSizeNumberLagrangianScalarFieldSource
:
    public LagrangianScalarFieldSource,
    private CloudLagrangianFieldSource<scalar>
{
public:

    // Public Enumerations

        //- Enumeration for the uniform size
        enum class uniformSize
        {
            number,
            surfaceArea,
            volume,
            mass
        };

        //- Names of the uniform size
        static const NamedEnum<uniformSize, 4> uniformSizeNames_;


private:

    // Private Data

        //- The size uniform to all parcels
        const uniformSize uniformSize_;


protected:

    // Protected Member Functions

        //- Calculate the necessary sizes
        void calcSizes
        (
            const LagrangianInjection& injection,
            const LagrangianSubMesh& subMesh,
            tmp<LagrangianSubScalarField>& size,
            const bool needv,
            tmp<LagrangianSubScalarField>& v,
            const bool needm,
            tmp<LagrangianSubScalarField>& m
        ) const;


public:

    //- Runtime type information
    TypeName("uniformSizeNumber");


    // Constructors

        //- Construct from internal field and dictionary
        uniformSizeNumberLagrangianScalarFieldSource
        (
            const regIOobject&,
            const dictionary& dict
        );

        //- Copy constructor setting the internal field reference
        uniformSizeNumberLagrangianScalarFieldSource
        (
            const uniformSizeNumberLagrangianScalarFieldSource&,
            const regIOobject&
        );

        //- Construct and return a clone setting the internal field reference
        virtual autoPtr<LagrangianScalarFieldSource> clone
        (
            const regIOobject& iIo
        ) const
        {
            return autoPtr<LagrangianScalarFieldSource>
            (
                new uniformSizeNumberLagrangianScalarFieldSource(*this, iIo)
            );
        }


    //- Destructor
    virtual ~uniformSizeNumberLagrangianScalarFieldSource();


    // Member Functions

        //- Return the sample size exponent
        label sampleQ() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
