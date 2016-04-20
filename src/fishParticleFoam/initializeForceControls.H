//******************************************************************************
//	Variable decleration for forceCalc.H.
//******************************************************************************


// Hard-code the patch name. Should be read from dictionary later...
const wordList motionPatches(1,"fish");
label forcePatchID 	= 	mesh.boundaryMesh().findPatchID("fish");
label inletPatchID 	= 	mesh.boundaryMesh().findPatchID("inlet");


// Create a field for viscous force.
Info << "Reading field viscousForce \n" << endl;
volVectorField viscousForce
(
	IOobject
	(
		"viscousForce",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh
);

// Create a field for total force.
Info << "Reading field totalForce \n" << endl;
volVectorField totalForce
(
	IOobject
	(
		"totalForce",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh
);


vectorField fishViscousForce;
vectorField fishTotalForce;
vectorField fishDragForce;

vectorField dragCheck;

scalar drag;
scalar inputPower;
scalar outputPower;
scalar froudeEfficiency;

