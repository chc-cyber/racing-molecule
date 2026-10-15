Clazz.declarePackage("java.net");
(function(){
var c$ = Clazz.declareType(java.net, "URLEncoder", null);
c$.encode = Clazz.defineMethod(c$, "encode", 
function(s){
return encodeURIComponent(s);
}, "~S");
})();
;//5.0.1-v7 Mon May 12 23:42:45 CDT 2025
