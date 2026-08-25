"use strict";

const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");

const applicationSource = fs.readFileSync(
  path.join(__dirname, "..", "magnetickp", "web", "static", "app.js"),
  "utf8",
);
const start = applicationSource.indexOf("function complex(");
const end = applicationSource.indexOf("function svgElement(");
assert.ok(start >= 0 && end > start, "band-math functions are present");
const bandMath = new Function(
  `const t = () => "invalid path"; ${applicationSource.slice(start, end)}; return { complex, evaluateTree, hermitianEigenvalues, parseBandPathValue, sampledBandPath };`,
)();
const { complex, evaluateTree, hermitianEigenvalues, parseBandPathValue, sampledBandPath } = bandMath;

function close(actual, expected, tolerance = 1e-8) {
  assert.ok(Math.abs(actual - expected) < tolerance, `${actual} != ${expected}`);
}

const diagonal = hermitianEigenvalues([
  [complex(1), complex(0)],
  [complex(0), complex(-1)],
]);
close(diagonal[0], -1);
close(diagonal[1], 1);

const complexHermitian = hermitianEigenvalues([
  [complex(0), complex(1, 1)],
  [complex(1, -1), complex(0)],
]);
close(complexHermitian[0], -Math.sqrt(2));
close(complexHermitian[1], Math.sqrt(2));

const treeValue = evaluateTree(
  {
    type: "add",
    arguments: [
      { type: "symbol", name: "x" },
      {
        type: "multiply",
        arguments: [
          { type: "complex", real: 0, imaginary: 1 },
          { type: "number", value: 2 },
        ],
      },
    ],
  },
  { x: 3 },
);
close(treeValue.real, 3);
close(treeValue.imaginary, 2);

const zeroSquared = evaluateTree(
  {
    type: "power",
    base: { type: "number", value: 0 },
    exponent: { type: "number", value: 2 },
  },
  {},
);
close(zeroSquared.real, 0);
close(zeroSquared.imaginary, 0);

const notebookPath = parseBandPathValue(
  "Γ:-2/3,1/3,0 -> K:0,0,0; K:0,0,0 -> M:0,1/2,0",
  3,
);
assert.equal(notebookPath.length, 2);
assert.deepEqual(notebookPath[0].start.coordinates, [-2 / 3, 1 / 3, 0]);
assert.deepEqual(notebookPath[1].end.coordinates, [0, 1 / 2, 0]);

const sampled = sampledBandPath(notebookPath, 2);
assert.equal(sampled.segments.length, 2);
assert.equal(sampled.segments[0].length, 3);
assert.equal(sampled.segments[1].length, 3);
assert.equal(sampled.totalDistance, 2);
assert.deepEqual(sampled.ticks.map((tick) => tick.label), ["Γ", "K", "M"]);

const disconnected = sampledBandPath(
  parseBandPathValue("A:0,0 -> B:1,0; C:0,1 -> D:1,1", 2),
  2,
);
assert.deepEqual(disconnected.ticks.map((tick) => tick.label), ["A", "B | C", "D"]);

console.log("MagneticKP browser band math: ok");
