Clazz.declarePackage("jme");
Clazz.load(["jme.JME"], "jme.JMEJmol", ["javax.swing.SwingUtilities"], function(){
var c$ = Clazz.decorateAsClass(function(){
this.parentWindow = null;
this.vwr = null;
this.headless = false;
Clazz.instantialize(this, arguments);}, jme, "JMEJmol", jme.JME, java.awt.event.WindowListener);
Clazz.defineMethod(c$, "setViewer", 
function(frame, vwr, parent, frameType){
this.parentWindow = parent;
this.vwr = vwr;
if (parent == null && frame == null && !"search".equals(frameType)) this.headless = vwr.headless;
if (!this.headless) {
if (frame == null) {
frame = getJmolFrame();
}this.setFrame(frame);
}this.initialize();
if (parent != null) {
if (vwr != null) vwr.getInchi(null, null, null);
javax.swing.SwingUtilities.invokeLater({
this.start();
});
}}, "javax.swing.JFrame,JV.Viewer,java.awt.Container,~S");
Clazz.overrideMethod(c$, "windowOpened", 
function(e){
}, "java.awt.event.WindowEvent");
Clazz.overrideMethod(c$, "windowClosing", 
function(e){
}, "java.awt.event.WindowEvent");
Clazz.overrideMethod(c$, "windowClosed", 
function(e){
if (this.myFrame != null) this.myFrame.setVisible(false);
}, "java.awt.event.WindowEvent");
Clazz.overrideMethod(c$, "windowIconified", 
function(e){
}, "java.awt.event.WindowEvent");
Clazz.overrideMethod(c$, "windowDeiconified", 
function(e){
}, "java.awt.event.WindowEvent");
Clazz.overrideMethod(c$, "windowActivated", 
function(e){
}, "java.awt.event.WindowEvent");
Clazz.overrideMethod(c$, "windowDeactivated", 
function(e){
}, "java.awt.event.WindowEvent");
Clazz.defineMethod(c$, "setFrameVisible", 
function(b){
if (this.myFrame != null) this.myFrame.setVisible(b);
}, "~B");
});
;//5.0.1-v7 Mon May 04 21:41:21 EDT 2026
