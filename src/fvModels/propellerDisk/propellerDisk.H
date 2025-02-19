/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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
    Foam::fv::propellerDisk

Description
    Disk momentum source which approximates a propeller based on a given
    propeller curve.

    Reference:
    \verbatim
        Hough, G. R., & Ordway, D. E. (1964).
        The generalized actuator disk.
        Developments in theoretical and applied mechanics, 2, 317-336.
    \endverbatim

Usage
    Example usage:
    \verbatim
    diskSource
    {
        type            propellerDisk;

        libs            ("libpropellerDisk.so");

        selectionMode   cellZone;
        cellZone        propeller;

        normal          (1 0 0);    // Normal direction of the propeller

        n               26.03;      // Rotation speed [1/s]

        dPropeller      0.203;      // Propeller diameter
        dHub            0.039179;   // Hub diameter

        propellerCurve
        {
            type table;

            //   J     Kt     Kq
            values
            (
                (0.10 (0.3267 0.03748))
                (0.15 (0.3112 0.03629))
                (0.20 (0.2949 0.03500))
                (0.25 (0.2777 0.03361))
                (0.30 (0.2598 0.03210))
                (0.35 (0.2410 0.03047))
                (0.40 (0.2214 0.02871))
                (0.45 (0.2010 0.02682))
                (0.50 (0.1798 0.02479))
                (0.55 (0.1577 0.02261))
                (0.60 (0.1349 0.02027))
                (0.65 (0.1112 0.01777))
                (0.70 (0.0867 0.01509))
                (0.75 (0.0614 0.01224))
                (0.80 (0.0353 0.00921))
            );
        }
    }
    \endverbatim

SourceFiles
    propellerDisk.C
    propellerDiskTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef propellerDisk_H
#define propellerDisk_H

#include "fvModel.H"
#include "logFile.H"
#include "fvCellSet.H"
#include "Function1.H"
#include "vector2D.H"
#include "forces.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class propellerDisk Declaration
\*---------------------------------------------------------------------------*/

class propellerDisk
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

        //- Propeller disk centre
        vector centre_;

        //- Propeller disk normal direction
        vector normal_;

        //- Rotation speed [1/s]
        scalar n_;

        //- Rotation direction (obtained from the sign of n_)
        scalar rotationDir_;

        //- Propeller diameter
        scalar dProp_;

        //- Hub diameter
        scalar dHub_;

        //- Propeller function
        autoPtr<Function1<vector2D>> propellerFunction_;

        mutable autoPtr<volVectorField::Internal> forcePtr_;

        //- Current force of the fluid on the propeller
        mutable vector force_;

        //- Current moment of the fluid on the propeller
        mutable vector moment_;

        //- Optional switch to enable logging of integral properties
        Switch log_;

        //- Optional log file
        mutable autoPtr<functionObjects::logFile> logFile_;


    // Protected Member Functions

        //- Read the model coefficients
        virtual void readCoeffs(const dictionary& dict);

        //- Return the propeller normal
        virtual vector normal() const
        {
            return normal_;
        }

        //- Return the rotation speed
        virtual scalar n() const
        {
            return n_;
        }

        //- No rotation speed correction
        virtual void correctn(const scalar T) const
        {}

        //- Computes the thickness of the disk in streamwise direction
        scalar diskThickness(const vector& centre) const;

        //- Return the normalised flow-rate through the disk
        scalar J(const vectorField& U, const vector& nHat) const;


private:

    // Private Member Functions

        //- Add resistance to the UEqn
        template<class AlphaFieldType, class RhoFieldType>
        void addActuationDiskAxialInertialResistance
        (
            vectorField& Usource,
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            const volVectorField::Internal& U
        ) const;


public:

    //- Runtime type information
    TypeName("propellerDisk");


    // Constructors

        //- Construct from components
        propellerDisk
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        propellerDisk(const propellerDisk&) = delete;


    //- Destructor
    virtual ~propellerDisk()
    {}


    // Member Functions

        //- Return the propeller centre
        virtual vector centre() const
        {
            return centre_;
        }

        //- Return the current force of the fluid on the propeller
        vector force() const
        {
            return force_;
        }

        //- Return the current moment of the fluid on the propeller
        vector moment() const
        {
            return moment_;
        }


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
        void operator=(const propellerDisk&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "propellerDiskTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
