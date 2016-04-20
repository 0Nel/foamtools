/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Application
    gen6DoF

Description
    Generate simple sinusoidal 6-DoF motion control-file.

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "vector.H"
#include "Vector2D.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "timeSelector.H"
#include "argList.H"
#include "IOobject.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
//#include "fvIOoptionList.H"



using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

	#include "addRegionOption.H"
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
	//#include "createNamedDynamicFvMesh.H"

	//const scalar rotAmplitude = runTime.controlDict().lookupOrDefault<scalar>("rotationAmplitude", 0);
	//const scalar rotFrequency = runTime.controlDict().lookupOrDefault<scalar>("rotationFrequency", 0);


/* keeping the code in case of need
	IOdictionary motionProperties
    (
        IOobject
        (
            "motionProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

	const scalar rotAmplitude (
    	readScalar(motionProperties.lookup("rotationAmplitude"))
	);

	const scalar rotFrequency (
		readScalar(motionProperties.lookup("rotationFrequency"))
	);
*/

	// read in the controlDict to get start time and end time
	IOdictionary controlDict
    (
        IOobject
        (
            "controlDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // End-time of the table
    const scalar endTime (readScalar(controlDict.lookup("endTime")));


	// read in the dynamicMeshDict
	IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

	// get the gen6DoF subDict
    dictionary gen6DoFDict (
    	dynamicMeshDict.subDict("generate6DoFDict")
    );

    const vector transAmplitude (gen6DoFDict.lookup("translationAmplitude"));

    const vector transFrequency (gen6DoFDict.lookup("translationFrequency"));

    const vector rotAmplitude (gen6DoFDict.lookup("rotationAmplitude"));

    const vector rotFrequency (gen6DoFDict.lookup("rotationFrequency"));

	Info << "translationAmplitude set to: "
		<< transAmplitude.x() << " " << transAmplitude.y() << " " << transAmplitude.z()
		<< endl;

	Info << "translationFrequency set to: "
		<< transFrequency.x() << " " << transFrequency.y() << " " << transFrequency.z()
		<< endl;

	Info << "rotationAmplitude set to: "
		<< rotAmplitude.x() << " " << rotAmplitude.y() << " " << rotAmplitude.z()
		<< endl;

	Info << "rotationFrequency set to: "
		<< rotFrequency.x() << " " << rotFrequency.y() << " " << rotFrequency.z()
		<< endl;


    // Number of entries in the table
    const scalar nTimes (readScalar(gen6DoFDict.lookup("nTimes")));

/*
    // Amplitude of the translation [m]
    const vector transAmp(0, 0, 0);

    // Frequency of the translation [rad/s]
    const vector transOmega(0.0, 0.0, 0.0);

    // Amplitude of the rotation [deg]
    const vector rotAmp(0, 45, 100);

    // Frequency of the rotation [rad/s]
    const vector rotOmega(1.0, 1.0, 1.0);
*/

	// define the output dictionary
    List<Tuple2<scalar, Vector2D<vector > > > timeValues(nTimes);

    const scalar timeStep = endTime / ( nTimes - 1);

    Info << "time step: " << timeStep <<endl;

    forAll(timeValues, i)
    {

        scalar t = (endTime*i)/(nTimes - 1);

		if(i != 0) {
	    	scalar tOld = (endTime*i)/(nTimes - 1) - (endTime / nTimes);
	   		Info << "current time: " << t << " old time: " << tOld << endl;
    	}

        timeValues[i].first() = t;

		// transalationsvektor, uninteressant
        timeValues[i].second()[0] = vector
        (
            transAmplitude.x()*Foam::sin(transFrequency.x()*t),
            transAmplitude.y()*Foam::sin(transFrequency.y()*t),
            transAmplitude.z()*Foam::sin(transFrequency.z()*t)
        );

        timeValues[i].second()[1] = vector
        (
            rotAmplitude.x()*Foam::sin(rotFrequency.x()*t),
            rotAmplitude.y()*Foam::sin(rotFrequency.y()*t),
            rotAmplitude.z()*Foam::sin(rotFrequency.z()*t)
        );
    }

    // write out the new dictionary:
    {
        OFstream dataFile("6DoF.dat");
        dataFile << timeValues << endl;
    }

// experimental stuff



	// Number of entries in the table
    //const scalar nTimes (readScalar(gen6DoFDict.lookup("nTimes")));

	scalar pi = 3.141592;

	// for testing purposes only
    vector a = vector (1, 0, 0);

	tensor rotX (
        1, 0, 					0,
        0, Foam::cos(pi),   -Foam::sin(pi),
        0, Foam::sin(pi),    Foam::cos(pi)
    );

	tensor rotY (
        Foam::cos(pi), 0, -Foam::sin(pi),
        0, 				  1,  0,
        Foam::sin(pi), 0,  Foam::cos(pi)
    );

	tensor rotZ (
        Foam::cos(pi), -Foam::sin(pi), 0,
        Foam::sin(pi),  Foam::cos(pi), 0,
        0, 				   0,                1
    );

	vector rotatedA = rotZ & a;

	Info << "a: " << a.x()
		<< " " << a.y()
		<< " " << a.z()
		<< endl
		<< "rotiert: "
		<< rotatedA.x()
		<< " " << rotatedA.y()
		<< " " << rotatedA.z()
		<< endl;

	List<Tuple2<scalar, Vector2D<vector > > > timeValuesExp(nTimes);

	forAll(timeValuesExp, i) {
		if (i == 0) { // if its the first time step, set everything to zero
			timeValuesExp[i].first() = 0;
			timeValuesExp[i].second()[0] = vector (0,0,0);
			timeValuesExp[i].second()[1] = vector (0,0,0);
			continue;
		}

		vector omegaOld = vector (
			timeValuesExp[i].second()[0].x(),
			timeValuesExp[i].second()[0].y(),
			timeValuesExp[i].second()[0].z()
		);


		vector omegaCur = vector (
			rotAmplitude.x() * Foam::sin(2 * pi * rotFrequency.x() * timeValuesExp[i].first()),
			rotAmplitude.y() * Foam::sin(2 * pi * rotFrequency.y() * timeValuesExp[i].first()),
			rotAmplitude.z() * Foam::sin(2 * pi * rotFrequency.z() * timeValuesExp[i].first())
		);

		// translation isn't of any interest (yet)
		timeValuesExp[i].second()[0] = vector (0,0,0);

		tensor rotXOld (
	        1, 0, 					0,
	        0, Foam::cos(pi),   -Foam::sin(pi),
	        0, Foam::sin(pi),    Foam::cos(pi)
	    );

		tensor rotYOld (
	        Foam::cos(pi), 0, -Foam::sin(pi),
	        0, 				  1,  0,
	        Foam::sin(pi), 0,  Foam::cos(pi)
	    );

		tensor rotZOld (
	        Foam::cos(pi), -Foam::sin(pi), 0,
	        Foam::sin(pi),  Foam::cos(pi), 0,
	        0, 				   0,                1
	    );

		tensor rotXCur (
	        1, 0, 					0,
	        0, Foam::cos(pi),   -Foam::sin(pi),
	        0, Foam::sin(pi),    Foam::cos(pi)
	    );

		tensor rotYCur (
	        Foam::cos(pi), 0, -Foam::sin(pi),
	        0, 				  1,  0,
	        Foam::sin(pi), 0,  Foam::cos(pi)
	    );

		tensor rotZCur (
	        Foam::cos(pi), -Foam::sin(pi), 0,
	        Foam::sin(pi),  Foam::cos(pi), 0,
	        0, 				   0,                1
	    );

	}

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
