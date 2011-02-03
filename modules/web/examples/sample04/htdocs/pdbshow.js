// pdbshow.js to be included in an html doc to toggle on and
//   off various representations fo proteins and ligands.
//
// These functions tied to onClick (or other) events.
// Toggle state is maintained by global variables here,
//  such as load_state, cartoon_state, etc.
//  These are implemented as a Hash object, defined in hash.js
//  which must be included in the html file that includes this pdbshow.js
//
var load_state = new Hash;
var zoom = ""; // let first molecule zoom
function get_protein(a) {
  // load protein, unless already loaded
  if (load_state[a] == undefined) {
    load = "load?filename=$PYMOL_PATH/modules/web/examples/data/" + a + ".pdb.gz" + zoom;
    cmd(load);
    zoom = "&zoom=0"; // subsequent loads shold not zoom
    cmd("orient?selection=bymol organic and elem n&animate=2");
    cmd("hide?representation=everything&selection=" + a);
    load_state[a] = true;
    cartoon_state[a] = false;
    surface_state[a] = false;
    nonbonded_state[a] = false;
    ribbon_state[a] = false;
    ligand_state[a] = false;
    active_site_state[a] = false;
  }
}
var cartoon_state = new Hash;
function toggle_cartoon(a) {
  // toggle cartoon display
  get_protein(a);
  var arg = "representation=cartoon&selection=" + a;
  if (cartoon_state[a] == true) {
    cmd("hide?" + arg);
    cartoon_state[a] = false;
  } else {
    cmd("show?" + arg);
    cartoon_state[a] = true;
  }
}
var surface_state = new Hash;
function toggle_surface(a) {
  // toggle surface display
  get_protein(a);
  var arg = "representation=surface&selection=polymer and (" + a + " within 5 of bymol (" + a + " and organic and elem n))";
  if (surface_state[a] == true) {
    cmd("hide?" + arg);
    surface_state[a] = false;
  } else {
    cmd("show?" + arg);
    cmd("set?name=two_sided_lighting");
    cmd("set?name=transparency&value=0.3");
    surface_state[a] = true;
  }
}
var ribbon_state = new Hash;
function toggle_ribbon(a) {
  // toggle polymer ribbon display
  get_protein(a);
  var arg = "representation=ribbon&selection=" + a + " and polymer";
  if (ribbon_state[a] == true) {
    cmd("hide?" + arg);
    ribbon_state[a] = false;
  } else {
    cmd("show?" + arg);
    ribbon_state[a] = true;
  }
}
var nonbonded_state = new Hash;
function toggle_nonbonded(a) {
  // toggle solvent/water display
  get_protein(a);
  var arg = "representation=nonbonded&selection=(" + a + " and solvent) within 5 of bymol (" + a + " and organic and elem n)";
  if (nonbonded_state[a] == true) {
    cmd("hide?" + arg);
    nonbonded_state[a] = false;
  } else {
    cmd("show?" + arg);
    nonbonded_state[a] = true;
  }
}
var ligand_state = new Hash;
function toggle_ligand(a) {
  // toggle organic/ligand display
  get_protein(a);
  var arg = "representation=sticks&selection=bymol (" + a + " and organic and elem n)";
  if (ligand_state[a] == true) {
    cmd("hide?" + arg);
    ligand_state[a] = false;
  } else {
    cmd("show?" + arg);
    ligand_state[a] = true;
  }
}
var active_site_state = new Hash;
function toggle_active_site(a) {
  // toggle active site/polymer near organic display
  get_protein(a);
  var arg = "representation=lines&selection=byres polymer and (" + a + " within 5 of bymol (" + a + " and organic and elem n))";
  if (active_site_state[a] == true) {
    cmd("hide?" + arg);
    active_site_state[a] = false;
  } else {
    cmd("show?" + arg);
    active_site_state[a] = true;
  }
}
var hb_state = new Hash;
function toggle_hb(a) {
  // toggle hydrogen bond display
  get_protein(a);
  var arg  = "name=" + a + "_hb";
  var arg2 = arg + "&selection1=" + a + " and polymer" + "&selection2=" + a + " and organic" + "&cutoff=3.2&mode=2";
  if (hb_state[a] == true) {
    cmd("delete?" + arg);
    hb_state[a] = false;
  } else {
    cmd("dist?" + arg2);
    hb_state[a] = true;
  }
}

function reinitialize() {
  // reinitialize PyMOL and my state variables
    //  if ( confirm("Really reinitialize") ) {
    cmd("reinitialize");
    var keys = load_state.keys();
    for (k in keys) {
      a = keys[k];
      load_state[a] = undefined;
      cartoon_state[a] = 0;
      surface_state[a] = 0;
      nonbonded_state[a] = 0;
      ribbon_state[a] = 0;
      ligand_state[a] = 0;
      active_site_state[a] = 0;
    }
    document.forms[0].reset();
    zoom = "";
    // }
}
