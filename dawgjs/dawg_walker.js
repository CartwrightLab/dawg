// Driver example code for dawgjs
var dawg = require('./dawg.js')

console.log("getRandomInt(0, 100)", dawg.getRandomInt(0, 100));
console.log("getRandomInt(0, 100)", dawg.getRandomInt(0, 100));
console.log("getRandomInt(0, 100)", dawg.getRandomInt(0, 100));

console.log("getRandomFloat(0, 100)", dawg.getRandomFloat(0, 100));
console.log("getRandomFloat(0, 100)", dawg.getRandomFloat(0, 100));
console.log("getRandomFloat(0, 100)", dawg.getRandomFloat(0, 100));

console.log("getDawgError('ahhhhh')", dawg.getDawgError('ahhhhh'));

var dawgWrapper = new dawg.DawgWrapper(4);
console.log("getRandom(0, 10000)", dawgWrapper.getRandom(0, 10000));
console.log("getRandom(0, 10000)", dawgWrapper.getRandom(0, 10000));
console.log("getRandom(0, 10000)", dawgWrapper.getRandom(0, 10000));
console.log("getRandom(0, 10000)", dawgWrapper.getRandom(0, 10000));
dawgWrapper.delete();
