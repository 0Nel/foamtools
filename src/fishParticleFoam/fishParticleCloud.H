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
    Foam::fishParticleCloud

Description
    A Cloud of solid particles

SourceFiles
    fishParticleCloudI.H
    fishParticleCloud.C
    fishParticleCloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef fishParticleCloud_H
#define fishParticleCloud_H

#include "Cloud.H"
#include "fishParticle.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class fishParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class fishParticleCloud
:
    public Cloud<fishParticle>
{
    // Private data

        const fvMesh& mesh_;

        IOdictionary particleProperties_;

        scalar rhop_;
        scalar e_;
        scalar mu_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fishParticleCloud(const fishParticleCloud&);

        //- Disallow default bitwise assignment
        void operator=(const fishParticleCloud&);


public:

    // Constructors

        //- Construct given mesh
        fishParticleCloud
        (
            const fvMesh&,
            const word& cloudName = "defaultCloud",
            bool readFields = true
        );


    // Member Functions

        // Access

            inline const fvMesh& mesh() const;

            inline scalar rhop() const;
            inline scalar e() const;
            inline scalar mu() const;


        // Edit

            //- Move the particles under the influence of the given
            //  gravitational acceleration
            void move(const dimensionedVector& g, int check);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fishParticleCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
