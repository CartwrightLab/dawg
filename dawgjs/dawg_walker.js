// Driver example code for dawgjs
var dawg = require('./dawg.js')

console.log("Hello dawg_walker here");

console.log("getRandomInt(0, 100)", dawg.getRandomInt(0, 100));
console.log("getRandomInt(0, 100)", dawg.getRandomInt(0, 100));
console.log("getRandomInt(0, 100)", dawg.getRandomInt(0, 100));

console.log("getRandomFloat(0, 100)", dawg.getRandomFloat(0, 100));
console.log("getRandomFloat(0, 100)", dawg.getRandomFloat(0, 100));
console.log("getRandomFloat(0, 100)", dawg.getRandomFloat(0, 100));

var dawgWrapper = new dawg.DawgWrapper(4444);
console.log("getDawgError('ahhhhh')", dawgWrapper.getDawgError('ahhhhh'));
console.log("getRandom(0, 10000)", dawgWrapper.getRandom(0, 10000));
console.log("getRandom(0, 10000)", dawgWrapper.getRandom(0, 10000));
console.log("getRandom(0, 10000)", dawgWrapper.getRandom(0, 10000));
console.log("getRandom(0, 10000)", dawgWrapper.getRandom(0, 10000));
console.log(dawgWrapper.getRnaSequence(60, "jc"));
dawgWrapper.delete();
