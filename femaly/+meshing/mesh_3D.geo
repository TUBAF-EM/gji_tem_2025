SetFactory("OpenCASCADE");

// Initialize points.
pt_list() = {};
ln_list() = {};
topo_pt_list() = {};

// Merge air-earth surface.
Include "mesh_3D_surf3D.geo";

// Get bb of surface for later identifictation of entities.
bb() = BoundingBox Surface {bnd_air_earth()};

/// Create domain, i.e. volume entity.
domain_r = 10000;
domain_c = {8.75, 0, 0};
shape = 1;
keep_air = 1;
earth_id = newv;
If (shape == 1)
	Sphere(earth_id) = {domain_c(0), domain_c(1), domain_c(2), domain_r, -Pi/2, Pi/2, 2*Pi};
Else
	Box(earth_id) = {domain_c(0)-domain_r, domain_c(1)-domain_r, domain_c(2)-domain_r, 2*domain_r, 2*domain_r, 2*domain_r};
EndIf

// Cut volume by spline surface (parts) and add points.
BooleanFragments{Volume{earth_id}; Delete;}{Surface{bnd_air_earth()}; Delete;}

// Catch errors.
tmp_vl() = Volume{:};
If (#tmp_vl() > 2)
	Error('Cutting volume by given surface results in more than the expected number of volumes.');
EndIf

// Set earth & air by BooleanFragments output numbering scheme.
If (shape == 1)
	air_id = tmp_vl(0);
	earth_id = tmp_vl(1);
Else
	air_id = tmp_vl(1);
	earth_id = tmp_vl(0);
EndIf

// Delete excess surfaces.
// Every needed surface is boundary of the volumes, hence delete every surface that isn't.
// Identify all surfaces of volumes.
domain_surfs() = Abs(Boundary{Volume{:};});
// Save all surfaces.
excess_surfs() = Surface{:};
// Leave only surfaces not part of volumes.
excess_surfs() -= domain_surfs();
Recursive Delete{Surface{excess_surfs()};}

// Add insitu entities.
pt_id = newp;
ln_id = newl;



// Identify remaining surfaces.
bnd_air() = Abs(Boundary{Volume{air_id};});
bnd_earth() = Abs(Boundary{Volume{earth_id};});
bnd_air_earth() = {};
For jj In {0:#bnd_earth()-1}
	For ii In {0:#bnd_air()-1}
		If (bnd_earth(jj) == bnd_air(ii))
			bnd_air_earth() += bnd_earth(jj);
		EndIf
	EndFor
EndFor
bnd_air() -= bnd_air_earth();
bnd_earth() -= bnd_air_earth();

// Remove or keep air domain.
If (keep_air == 1)
	Physical Volume("air", 1) = air_id;
	Physical Surface("bnd_air", 1) = bnd_air();
Else
	Recursive Delete{Volume{air_id};}
EndIf

// Set physical entities.
Physical Volume("halfspace", 2) = earth_id;
Physical Surface("bnd_earth", 2) = bnd_earth();
Physical Surface("bnd_air_earth", 3) = bnd_air_earth();

/// Meshing.
// Initialize meshing parameter.
size_at_point = 0.5;
size_at_wire = 1.25;
pad = 0.7;

// Cell sizes at points.
Field[1] = Distance;
Field[1].NodesList = {pt_list()};
Field[10] = Threshold;
Field[10].IField = 1;
Field[10].LcMin = size_at_point;
Field[10].DistMin = 2*size_at_point;
Field[10].LcMax = pad*domain_r;
Field[10].DistMax = domain_r;
MeshSize{pt_list()} = size_at_point; // FIXME: why additionaly required?

// Cell sizes at lines.
Field[2] = Distance;
Field[2].EdgesList = {ln_list()};
Field[20] = Threshold;
Field[20].IField = 2;
Field[20].LcMin = size_at_wire;
Field[20].DistMin = size_at_wire;
Field[20].LcMax = pad*domain_r;
Field[20].DistMax = domain_r;

// Take the min of all constraints.
Field[100] = Min;
Field[100].FieldsList = {10, 20};
Background Field = {100};
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.SaveParametric = 0;

/// User geo code.
Field[6] = Box;
Field[6].VIn = 1.5;
Field[6].VOut = domain_r;
Field[6].XMin = 0;
Field[6].XMax = 20;
Field[6].YMin = -5;
Field[6].YMax = 5;
Field[6].ZMin = -5;
Field[6].ZMax = 0;
Field[6].Thickness = 10;
Field[200] = Min;
Field[200].FieldsList = {10, 20, 6};
Background Field = {200};


// Run Gmsh.
// FIXME: during optimization seg. faults may occur
// Mesh.OptimizeNetgen = 2;
Mesh 3;

