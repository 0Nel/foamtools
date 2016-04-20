/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::fishParticle

Description
    Simple solid spherical particle class with one-way coupling with the
    continuous phase.

SourceFiles
    fishParticleI.H
    fishParticle.C
    fishParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef fishParticle_H
#define fishParticle_H

#include "Particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fishParticleCloud;

/*---------------------------------------------------------------------------*\
                           Class fishParticle Declaration
\*---------------------------------------------------------------------------*/

class fishParticle
:
    public Particle<fishParticle>
{
    // Private member data

        //- Diameter
        scalar d_;

        //- Velocity of parcel
        vector U_;


public:

    friend class Cloud<fishParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackData
    {
        //- Reference to the cloud containing this particle
        fishParticleCloud& spc_;

        // Interpolators for continuous phase fields

            const interpolationCellPoint<scalar>& rhoInterp_;
            const interpolationCellPoint<vector>& UInterp_;
            const interpolationCellPoint<scalar>& nuInterp_;

        //- Local gravitational or other body-force acceleration
        const vector& g_;


    public:

        bool switchProcessor;
        bool keepParticle;


        // Constructors

            inline trackData
            (
                fishParticleCloud& spc,
                const interpolationCellPoint<scalar>& rhoInterp,
                const interpolationCellPoint<vector>& UInterp,
                const interpolationCellPoint<scalar>& nuInterp,
                const vector& g
            );


        // Member functions

            inline fishParticleCloud& spc();

            inline const interpolationCellPoint<scalar>& rhoInterp() const;

            inline const interpolationCellPoint<vector>& UInterp() const;

            inline const interpolationCellPoint<scalar>& nuInterp() const;

            inline const vector& g() const;
    };


    // Constructors

        //- Construct from components
        inline fishParticle
        (
            const Cloud<fishParticle>& c,
            const vector& position,
            const label celli,
            const scalar m,
            const vector& U
        );

        //- Construct from Istream
        fishParticle
        (
            const Cloud<fishParticle>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<fishParticle> clone() const
        {
            return autoPtr<fishParticle>(new fishParticle(*this));
        }


    // Member Functions

        // Access

            //- Return diameter
            inline scalar d() const;

            //- Return velocity
            inline const vector& U() const;

            //- The nearest distance to a wall that
            //  the particle can be in the n direction
            inline scalar wallImpactDistance(const vector& n) const;


        // Tracking

            //- Move
            bool move(trackData&);


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch
            (
                const polyPatch&,
                fishParticle::trackData& td,
                const label patchI
            );

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions without trackData
            bool hitPatch
            (
                const polyPatch& p,
                int& td,
                const label patchI
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                fishParticle::trackData& td
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch without trackData
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch
            (
                const wallPolyPatch&,
                fishParticle::trackData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            //- without trackData
            void hitWallPatch
            (
                const wallPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a polyPatch
            void hitPatch
            (
                const polyPatch&,
                fishParticle::trackData& td
            );

            //- Overridable function to handle the particle hitting a polyPatch
            //- without trackData
            void hitPatch
            (
                const polyPatch&,
                int&
            );

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            void transformProperties
            (
                const tensor& T
            );

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            void transformProperties
            (
                const vector& separation
            );


    // I-O

        static void readFields(Cloud<fishParticle>& c);

        static void writeFields(const Cloud<fishParticle>& c);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const fishParticle&);
};


template<>
inline bool contiguous<fishParticle>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fishParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
