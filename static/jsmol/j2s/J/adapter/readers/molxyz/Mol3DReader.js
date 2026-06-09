Clazz.declarePackage("J.adapter.readers.molxyz");
Clazz.load(["J.adapter.readers.molxyz.MolReader"], "J.adapter.readers.molxyz.Mol3DReader", null, function(){
var c$ = Clazz.declareType(J.adapter.readers.molxyz, "Mol3DReader", J.adapter.readers.molxyz.MolReader);
Clazz.overrideMethod(c$, "initializeReader", 
function(){
this.allow2D = false;
});
});
;//5.0.1-v7 Wed Jun 03 09:22:38 CDT 2026
