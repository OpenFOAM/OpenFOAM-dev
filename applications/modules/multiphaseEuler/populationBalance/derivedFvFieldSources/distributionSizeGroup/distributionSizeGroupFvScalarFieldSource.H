/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::distributionSizeGroupFvScalarFieldSource

Description
    Source condition for size-group fractions that sets values corresponding
    to a given distribution

Usage
    Example of the source condition specification:
    \verbatim
    <sourceName>
    {
        type        distributionSizeGroup;
        distribution
        {
            type        RosinRammler;
            Q           0;
            min         1 [mm];
            max         15 [mm];
            d           14 [mm];
            n           2 [];
        }
    }
    \endverbatim

SourceFiles
    distributionSizeGroupFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef distributionSizeGroupFvScalarFieldSource_H
#define distributionSizeGroupFvScalarFieldSource_H

#include "fvFieldSources.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class distributionSizeGroupFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class distributionSizeGroupFvScalarFieldSource
:
    public fvScalarFieldSource
{
private:

    // Private Data

        //- Distribution
        const autoPtr<distribution> distribution_;

        //- Cached allocation coefficient value. Calculated on demand.
        mutable autoPtr<scalar> etaPtr_;


public:

    //- Runtime type information
    TypeName("distributionSizeGroup");


    // Constructors

        //- Construct from internal field and dictionary
        distributionSizeGroupFvScalarFieldSource
        (
            const DimensionedField<scalar, volMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        distributionSizeGroupFvScalarFieldSource
        (
            const distributionSizeGroupFvScalarFieldSource&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<fvScalarFieldSource> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return autoPtr<fvScalarFieldSource>
            (
                new distributionSizeGroupFvScalarFieldSource
                (
                    *this,
                    iF
                )
            );
        }


    //- Destructor
    virtual ~distributionSizeGroupFvScalarFieldSource();


    // Member Functions

        //- Return the source value
        virtual tmp<scalarField> sourceValue(const fvSource&) const;

        //- Return the internal coefficient
        virtual tmp<scalarField> internalCoeff(const fvSource&) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
