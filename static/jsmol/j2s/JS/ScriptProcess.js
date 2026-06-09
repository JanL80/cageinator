Clazz.declarePackage("JS");
(function(){
var c$ = Clazz.decorateAsClass(function(){
this.processName = null;
this.context = null;
Clazz.instantialize(this, arguments);}, JS, "ScriptProcess", null);
Clazz.makeConstructor(c$, 
function(name, context){
this.processName = name;
this.context = context;
}, "~S,JS.ScriptContext");
})();
;//5.0.1-v7 Wed Jun 03 09:22:38 CDT 2026
