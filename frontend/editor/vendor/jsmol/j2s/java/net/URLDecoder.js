Clazz.declarePackage("java.net");
(function(){
var c$ = Clazz.declareType(java.net, "URLDecoder", null);
c$.decode = Clazz.defineMethod(c$, "decode", 
function(s){
return decodeURIComponent(s);
}, "~S");
})();
;//5.0.1-v7 Mon May 12 23:42:45 CDT 2025
