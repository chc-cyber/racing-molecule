Clazz.declarePackage("J.adapter.writers");
Clazz.load(["J.adapter.writers.CIFWriter"], "J.adapter.writers.SCIFWriter", null, function(){
var c$ = Clazz.decorateAsClass(function(){
this.scifInfo = null;
this.headerKeys = "\n_space_group_spin.transform_spinframe_P_abc SSG_SF\n_space_group_spin.collinear_direction SSG_CPD\n_space_group_spin.coplanar_perp_uvw SSG_CPU\n_space_group_spin.rotation_axis_xyz SSG_ROTAXIS\n_space_group_spin.rotation_angle SSG_ROTANGE\n_space_group_spin.number_SpSG_Chen  SSG_NUMBER\n_space_group_spin.name_SpSG_Chen    SSG_NAME\n\n";
this.operationKeys = "\nloop_\n_space_group_symop_spin_operation.id\n_space_group_symop_spin_operation.xyzt\n_space_group_symop_spin_operation.uvw_id\n";
this.latticeKeys = "\nloop_\n_space_group_symop_spin_lattice.id\n_space_group_symop_spin_lattice.xyzt\n_space_group_symop_spin_lattice.uvw_id\n";
this.upartKeys = "\nloop_\n_space_group_symop_spin_Upart.id\n_space_group_symop_spin_Upart.uvw\n";
this.momentKeys = "\nloop_\n_atom_site_spin_moment.label\n_atom_site_spin_moment.axis_u\n_atom_site_spin_moment.axis_v\n_atom_site_spin_moment.axis_w\n_atom_site_spin_moment.symmform_uvw\n_atom_site_spin_moment.magnitude\n";
Clazz.instantialize(this, arguments);}, J.adapter.writers, "SCIFWriter", J.adapter.writers.CIFWriter);
Clazz.makeConstructor(c$, 
function(){
Clazz.superConstructor(this, J.adapter.writers.SCIFWriter);
this.isCIF2 = true;
});
Clazz.defineMethod(c$, "prepareAtomSet", 
function(bs){
bs = Clazz.superCall(this, J.adapter.writers.SCIFWriter, "prepareAtomSet", [bs]);
this.scifInfo = this.vwr.getModelSetAuxiliaryInfoForAtoms(bs).get("scifInfo");
return bs;
}, "JU.BS");
Clazz.overrideMethod(c$, "writeHeader", 
function(sb){
}, "JU.SB");
Clazz.overrideMethod(c$, "writeOperations", 
function(sb){
this.scifInfo.get("");
sb.append(this.operationKeys);
sb.append(this.latticeKeys);
sb.append(this.upartKeys);
}, "JU.SB");
Clazz.defineMethod(c$, "writeAtomSite", 
function(sb){
sb.append(this.momentKeys);
for (var i = this.bsOut.nextSetBit(0), p = 1; i >= 0; i = this.bsOut.nextSetBit(i + 1)) {
var a = this.atoms[i];
var m = a.getVibrationVector();
if (m.magMoment == 0) continue;
sb.appendI(p++);
this.append3(sb, m);
sb.append(" ?");
sb.append(" ").appendF(m.length());
sb.append("\n");
}
return Clazz.superCall(this, J.adapter.writers.SCIFWriter, "writeAtomSite", [sb]);
}, "JU.SB");
});
;//5.0.1-v7 Mon Jun 09 19:31:20 CEST 2025
