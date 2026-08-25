"use strict";

const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");

const applicationSource = fs.readFileSync(
  path.join(__dirname, "..", "magnetickp", "web", "static", "app.js"),
  "utf8",
);
const start = applicationSource.indexOf("function unicodeCrystallographicSymbol(");
const end = applicationSource.indexOf("function selectedCatalogPoint(");
assert.ok(start >= 0 && end > start, "catalog symbol formatter is present");
const unicodeCrystallographicSymbol = new Function(
  `${applicationSource.slice(start, end)}; return unicodeCrystallographicSymbol;`,
)();

assert.equal(unicodeCrystallographicSymbol("P4_12_12"), "P4₁2₁2");
assert.equal(unicodeCrystallographicSymbol("P6_3/mmc"), "P6₃/mmc");
assert.equal(unicodeCrystallographicSymbol("I4_1/acd"), "I4₁/acd");
assert.equal(unicodeCrystallographicSymbol("P2_{1}/c"), "P2₁/c");
assert.equal(unicodeCrystallographicSymbol("P\\bar{4}2_1m"), "P4̅2₁m");
assert.equal(unicodeCrystallographicSymbol("Pm\\bar{3}m"), "Pm3̅m");
assert.equal(unicodeCrystallographicSymbol("P4_{1}2_{1}21'"), "P4₁2₁21′");
assert.equal(unicodeCrystallographicSymbol("P_{C}4_{1}2_{1}2"), "P₍C₎4₁2₁2");
assert.equal(unicodeCrystallographicSymbol("P_{s}\\bar{1}"), "P₍s₎1̅");
assert.equal(unicodeCrystallographicSymbol("", "P1"), "P1");

console.log("MagneticKP catalog symbols: ok");
