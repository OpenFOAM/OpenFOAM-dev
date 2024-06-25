/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::fv::actuationDisk

Description
    Actuation disk source

    Constant values for momentum source for actuation disk
    \f[
        T = 2 \rho A (\hat{d}\dot U_{o})^2 a (1-a) \hat{d}
    \f]

    where:
    \vartable
        A    | Disk area
        dHat | Unit disk direction
        U_o  | Upstream velocity
        a    | 1 - Cp/Ct
        Cp   | Power coefficient
        Ct   | Thrust coefficient
    \endvartable

Usage
    Example usage:
    \verbatim
    actuationDisk1
    {
        type            actuationDisk;

        select          cellSet;
        cellSet         actuationDisk1;

        diskDir         (-1 0 0);   // Disk direction
        Cp              0.1;        // Power coefficient
        Ct              0.5;        // Thrust coefficient
        diskArea        5.0;        // Disk area
        upstreamPoint   (0 0 0);    // Upstream point
    }
    \endverbatim

SourceFiles
    actuationDisk.C
    actuationDiskTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef actuationDisk_H
#define actuationDisk_H

#include "fvModel.H"
#include "fvCellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class actuationDisk Declaration
\*---------------------------------------------------------------------------*/

class actuationDisk
:
    public fvModel
{
protected:

    // Protected Data

        //- The set of cells the fvConstraint applies to
        fvCellSet set_;

        //- The name of the phase to which this fvModel applies
        word phaseName_;

        //- Name of the velocity field
        word UName_;

        //- Disk area normal
        vector diskDir_;

        //- Power coefficient
        scalar Cp_;

        //- Thrust coefficient
        scalar Ct_;

        //- Disk area
        scalar diskArea_;

        //- Upstream point sample
        point upstreamPoint_;

        //- Upstream cell ID
        label upstreamCellId_;


private:

    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- Add resistance to the UEqn
        template<class AlphaFieldType, class RhoFieldType>
        inline void addActuationDiskAxialInertialResistance
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;


public:

    //- Runtime type information
    TypeName("actuationDisk");


    // Constructors

        //- Construct from components
        actuationDisk
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        actuationDisk(const actuationDisk&) = delete;


    //- Destructor
    virtual ~actuationDisk()
    {}


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


         // Add explicit and implicit contributions

            //- Source term to momentum equation
            virtual void addSup
            (
                const volVectorField& U,
                fvMatrix<vector>& eqn
            ) const;

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volVectorField& U,
                fvMatrix<vector>& eqn
            ) const;

            //- Explicit and implicit sources for phase equations
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volVectorField& U,
                fvMatrix<vector>& eqn
            ) const;


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const actuationDisk&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
