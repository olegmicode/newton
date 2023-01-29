import Newton from "./libs/newton.js";

const fun = (x) => {
  return [
    x[0] + 0.5 * Math.pow(x[0] - x[1], 3) - 1.0,
    0.5 * Math.pow(x[1] - x[0], 3) + x[1],
  ];
};

const sol = Newton.Solve(fun, [0, 0]);
console.log(sol); // [0.84116396, 0.15883641] expected result
